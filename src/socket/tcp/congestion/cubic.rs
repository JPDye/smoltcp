use crate::{socket::tcp::RttEstimator, time::Instant};

use super::Controller;

// Constants for the Cubic congestion control algorithm.
// See RFC 8312.
const BETA_CUBIC: f64 = 0.7;
const C: f64 = 0.4;

const DEFAULT_MSS: usize = 1024;

#[derive(Debug)]
#[cfg_attr(feature = "defmt", derive(defmt::Format))]
pub struct Cubic {
    w_max: usize, // window size prior to loss
    cwnd: usize,
    min_cwnd: usize,
    ssthresh: usize,
    rwnd: usize,

    recovery_start: Option<Instant>,
    in_fast_recovery: bool,
}

impl Cubic {
    pub fn new() -> Cubic {
        Cubic {
            w_max: DEFAULT_MSS * 2,
            cwnd: DEFAULT_MSS * 2,
            min_cwnd: DEFAULT_MSS * 2,
            rwnd: 64 * DEFAULT_MSS,
            ssthresh: usize::MAX,

            recovery_start: None,
            in_fast_recovery: false,
        }
    }
}

impl Controller for Cubic {
    fn window(&self) -> usize {
        self.cwnd
    }

    fn on_ack(&mut self, now: Instant, len: usize, _in_flight: usize, _rtt: &RttEstimator) {
        // First new-data-ack exits fast recovery and deflates `cwnd`
        if self.in_fast_recovery {
            self.in_fast_recovery = false;
            self.cwnd = self.ssthresh;
            return;
        } else if self.cwnd < self.ssthresh {
            // Slow start: increase `cwnd` by 1 MSS per ACK.
            self.cwnd = self
                .cwnd
                .saturating_add(len.min(self.min_cwnd))
                .min(self.rwnd)
                .max(self.min_cwnd);
            return;
        }

        // Congestion avoidance: use cubic profile to advance window
        let recovery_start = self
            .recovery_start
            .expect("can't enter CA without having experienced loss");

        // Elapsed time since the start of the recovery phase.
        let t = now.total_millis() - recovery_start.total_millis();
        if t < 0 {
            return;
        }

        // RFC defines C in segments/sec^3; so scale to bytes/sec^3 since w_max and cwnd are bytes.
        let c_bytes = C * self.min_cwnd as f64;

        // K = (w_max * (1 - beta) / C)^(1/3)
        let k3 = ((self.w_max as f64) * (1.0 - BETA_CUBIC)) / c_bytes;
        let k = if let Some(k) = cube_root(k3) {
            k
        } else {
            return;
        };

        // cwnd = C(T - K)^3 + w_max
        let s = t as f64 / 1000.0 - k;
        let s = s * s * s;
        let cwnd = c_bytes * s + self.w_max as f64;

        self.cwnd = (cwnd as usize).max(self.min_cwnd).min(self.rwnd);
    }

    fn on_dup_ack(&mut self, _now: Instant, len: usize, _in_flight: usize) {
        if self.in_fast_recovery {
            self.cwnd = self
                .cwnd
                .saturating_add(len)
                .min(self.rwnd)
                .max(self.min_cwnd);
        }
    }

    fn on_loss(&mut self, now: Instant, _in_flight: usize) {
        // Only cut window size on first entrance to fast recovery.
        if !self.in_fast_recovery {
            // TODO: Make this optional?
            // RFC recommends (SHOULD) disabling if only a single CUBIC flow is on a network.
            //
            // RFC 9483.4.7: Fast Convergence
            // If loss happened at a smaller cwnd than before, it indicates a new flow.
            // Reduce the cubic plateau more than usual to create headroom.
            self.w_max = if self.cwnd < self.w_max {
                ((self.cwnd as f64) * (1.0 + BETA_CUBIC) / 2.0) as usize
            } else {
                self.cwnd
            };

            self.ssthresh = (((self.cwnd as f64) * BETA_CUBIC) as usize).max(2 * self.min_cwnd);
            self.cwnd = self
                .ssthresh
                .saturating_add(3 * self.min_cwnd)
                .min(self.rwnd);

            self.recovery_start = Some(now);
            self.in_fast_recovery = true;
        }
    }

    fn on_rto(&mut self, now: Instant, in_flight: usize) {
        self.w_max = self.cwnd;
        self.ssthresh = (in_flight >> 1).max(2 * self.min_cwnd);
        self.cwnd = self.min_cwnd;

        self.recovery_start = Some(now);
        self.in_fast_recovery = false
    }

    fn set_mss(&mut self, mss: usize) {
        self.min_cwnd = mss;
    }

    fn set_remote_window(&mut self, remote_window: usize) {
        if self.rwnd < remote_window {
            self.rwnd = remote_window;
        }
    }
}

#[inline]
fn abs(a: f64) -> f64 {
    if a < 0.0 { -a } else { a }
}

/// Calculate cube root by using the Newton-Raphson method.
fn cube_root(a: f64) -> Option<f64> {
    if a <= 0.0 {
        return None;
    }

    let (tolerance, init) = if a < 1_000.0 {
        (1.0, 8.879040017426005) // cube_root(700.0)
    } else if a < 1_000_000.0 {
        (5.0, 88.79040017426004) // cube_root(700_000.0)
    } else if a < 1_000_000_000.0 {
        (50.0, 887.9040017426004) // cube_root(700_000_000.0)
    } else if a < 1_000_000_000_000.0 {
        (500.0, 8879.040017426003) // cube_root(700_000_000_000.0)
    } else if a < 1_000_000_000_000_000.0 {
        (5000.0, 88790.40017426001) // cube_root(700_000_000_000.0)
    } else {
        (50000.0, 887904.0017426) // cube_root(700_000_000_000_000.0)
    };

    let mut x = init; // initial value
    let mut n = 20; // The maximum iteration
    loop {
        let next_x = (2.0 * x + a / (x * x)) / 3.0;
        if abs(next_x - x) < tolerance {
            return Some(next_x);
        }
        x = next_x;

        if n == 0 {
            return Some(next_x);
        }

        n -= 1;
    }
}

#[cfg(test)]
mod test {
    use crate::{socket::tcp::RttEstimator, time::Instant};

    use super::*;

    const MSS: usize = 1024;

    fn ack(cubic: &mut Cubic, len: usize, now: Instant) {
        cubic.on_ack(now, len, cubic.window().saturating_sub(MSS), &rtte())
    }

    fn rtte() -> RttEstimator {
        RttEstimator::default()
    }

    #[test]
    fn congestion_avoidance_works() {
        let mut cubic = Cubic::new();
        cubic.set_mss(MSS);
        cubic.w_max = MSS * 32;

        // Post-fast-recovery state: cwnd = ssthresh ≈ w_max * beta.
        cubic.cwnd = (MSS * 32 * 7) / 10;
        cubic.ssthresh = cubic.cwnd;
        cubic.recovery_start = Some(Instant::from_millis(0));

        // CA at small time intervals should grow by less than 1 MSS per ACK.
        for i in 1..10 {
            let initial_cwnd = cubic.window();
            ack(&mut cubic, MSS, Instant::from_millis(i));
            assert!(cubic.window() < initial_cwnd + MSS);
        }

        // CA approaches w_max as t approaches K, and exceeds it past K.
        let pre = cubic.window();
        ack(&mut cubic, MSS, Instant::from_millis(3000));
        assert!(cubic.window() >= cubic.w_max);
        assert!(cubic.window() > pre);

        // CA should cap at the receive window.
        ack(&mut cubic, MSS, Instant::from_millis(100_000));
        assert_eq!(cubic.window(), cubic.rwnd);
    }

    #[test]
    fn fast_recovery_works() {
        let mut cubic = Cubic::new();
        cubic.set_mss(MSS);
        cubic.cwnd = MSS * 32;

        // duplicate ACKs before fast recovery should do nothing
        let initial_cwnd = cubic.window();
        for _ in 0..3 {
            cubic.on_dup_ack(Instant::from_millis(0), MSS, initial_cwnd);
        }
        assert_eq!(cubic.window(), initial_cwnd);

        // we enter fast recovery upon minor loss (three duplicate ACKs).
        // ssthresh = cwnd * beta_cubic, cwnd = ssthresh + 3*MSS, recovery_start = now.
        // w_max = cwnd since the prior w_max (initial 2*MSS) is below cwnd.
        let expected_ssthresh = (initial_cwnd as f64 * BETA_CUBIC) as usize;
        cubic.on_loss(Instant::from_millis(0), initial_cwnd / 2);
        assert_eq!(cubic.ssthresh, expected_ssthresh);
        assert_eq!(cubic.cwnd, expected_ssthresh + 3 * MSS);
        assert_eq!(cubic.w_max, initial_cwnd);
        assert!(cubic.in_fast_recovery);
        assert_eq!(cubic.recovery_start, Some(Instant::from_millis(0)));

        // in fast recovery, each dup-ACK should increase the cwnd by 1 MSS
        let initial_cwnd = cubic.window();
        for i in 0..3 {
            for _ in 0..3 {
                let initial_cwnd = cubic.window();
                cubic.on_dup_ack(Instant::from_millis(i), MSS, initial_cwnd);
                assert_eq!(cubic.window(), initial_cwnd + MSS);
            }

            // multiple loss events (trip-dup-ack) should not trigger additional fast recovery reductions
            let initial_cwnd = cubic.window();
            let initial_ssthresh = cubic.ssthresh;
            let initial_w_max = cubic.w_max;
            cubic.on_loss(Instant::from_millis(i), initial_cwnd);
            assert_eq!(cubic.window(), initial_cwnd);
            assert_eq!(cubic.ssthresh, initial_ssthresh);
            assert_eq!(cubic.w_max, initial_w_max);
        }
        assert_eq!(cubic.window(), initial_cwnd + MSS * 9);

        // a non-duplicate ACK exits fast recovery and deflates cwnd to ssthresh
        ack(&mut cubic, MSS, Instant::from_millis(10));
        assert_eq!(cubic.window(), cubic.ssthresh);
        assert!(!cubic.in_fast_recovery);
    }

    #[test]
    fn slow_start_works() {
        let mut cubic = Cubic::new();
        cubic.set_mss(MSS);
        cubic.cwnd = MSS * 32;
        cubic.ssthresh = MSS * 16;

        // we enter slow start upon major loss (an RTO)
        // window resets to MSS, ssthresh becomes half the in-flight bytes,
        // and recovery_start is updated so any later CA uses a fresh epoch.
        let initial_cwnd = cubic.window();
        let inflight = initial_cwnd;
        cubic.on_rto(Instant::from_millis(0), inflight);
        assert_eq!(cubic.ssthresh, inflight / 2);
        assert_eq!(cubic.window(), MSS);
        assert!(!cubic.in_fast_recovery);
        assert_eq!(cubic.recovery_start, Some(Instant::from_millis(0)));
        assert_eq!(cubic.w_max, initial_cwnd);

        // slow start grows by at most the MSS per ack
        let initial_cwnd = cubic.window();
        for i in 0..10 {
            let initial_cwnd = cubic.window();
            let now = Instant::from_millis(i);
            ack(&mut cubic, MSS * 2, now);
            assert_eq!(cubic.window(), initial_cwnd + MSS);
        }
        assert_eq!(cubic.window(), initial_cwnd + MSS * 10);

        // slow start uses the number of ACKed bytes if they're less than the MSS
        let initial_cwnd = cubic.window();
        for i in 0..10 {
            let initial_cwnd = cubic.window();
            let now = Instant::from_millis(10 + i);
            ack(&mut cubic, MSS / 2, now);
            assert_eq!(cubic.window(), initial_cwnd + MSS / 2);
        }
        assert_eq!(cubic.window(), initial_cwnd + MSS / 2 * 10);

        // slow start transitions to congestion avoidance at ssthresh
        let initial_cwnd = cubic.window();
        cubic.ssthresh = initial_cwnd + MSS;
        ack(&mut cubic, MSS, Instant::from_millis(30));
        assert_eq!(cubic.window(), initial_cwnd + MSS);
        assert_eq!(cubic.ssthresh, initial_cwnd + MSS);
    }

    #[test]
    fn progress_to_ca_via_rto() {
        let mut cubic = Cubic::new();
        cubic.set_mss(MSS);

        let mut time = 0;

        // slow start from default state
        let initial_cwnd = cubic.window();
        for _ in 0..30 {
            time += 1;
            ack(&mut cubic, MSS, Instant::from_millis(time));
        }
        assert_eq!(cubic.window(), initial_cwnd + MSS * 30);
        assert!(cubic.window() < cubic.ssthresh);

        // rto: cwnd resets to MSS, ssthresh becomes half in-flight bytes
        let rto_cwnd = cubic.window();
        cubic.on_rto(Instant::from_millis(time), rto_cwnd);
        assert_eq!(cubic.window(), MSS);
        assert_eq!(cubic.ssthresh, rto_cwnd / 2);
        assert_eq!(cubic.w_max, rto_cwnd);

        // slow start again until cwnd reaches new ssthresh
        while cubic.window() < cubic.ssthresh {
            time += 1;
            let initial_cwnd = cubic.window();
            ack(&mut cubic, MSS, Instant::from_millis(time));
            assert_eq!(cubic.window(), initial_cwnd + MSS);
        }
        assert_eq!(cubic.window(), cubic.ssthresh);

        // ca: subsequent ACKs follow the cubic curve
        time += 1;
        let initial_cwnd = cubic.window();
        ack(&mut cubic, MSS, Instant::from_millis(time));
        assert!(cubic.window() > initial_cwnd);
    }

    #[test]
    fn progress_to_ca_via_loss() {
        let mut cubic = Cubic::new();
        cubic.set_mss(MSS);

        let mut time = 0;

        // slow start from default state
        let initial_cwnd = cubic.window();
        for _ in 0..30 {
            time += 1;
            ack(&mut cubic, MSS, Instant::from_millis(time));
        }
        assert_eq!(cubic.window(), initial_cwnd + MSS * 30);
        assert!(cubic.window() < cubic.ssthresh);

        // dup ACKs: ssthresh = cwnd * beta, cwnd = ssthresh + 3*MSS, recovery_start = now
        time += 1;
        let loss_cwnd = cubic.window();
        let expected_ssthresh = (loss_cwnd as f64 * BETA_CUBIC) as usize;
        cubic.on_loss(Instant::from_millis(time), loss_cwnd);
        assert_eq!(cubic.ssthresh, expected_ssthresh);
        assert_eq!(cubic.window(), expected_ssthresh + 3 * MSS);
        assert!(cubic.in_fast_recovery);
        assert_eq!(cubic.recovery_start, Some(Instant::from_millis(time)));

        // inflate cwnd on each duplicate ACK
        for _ in 0..9 {
            time += 1;
            let initial_cwnd = cubic.window();
            cubic.on_dup_ack(Instant::from_millis(time), MSS, cubic.cwnd);
            assert_eq!(cubic.window(), initial_cwnd + MSS);
        }

        // non-duplicate ACK deflates cwnd to ssthresh
        time += 1;
        ack(&mut cubic, MSS, Instant::from_millis(time));
        assert_eq!(cubic.window(), expected_ssthresh);
        assert!(!cubic.in_fast_recovery);

        // ca: subsequent ACKs follow the cubic curve
        time += 1;
        let initial_cwnd = cubic.window();
        ack(&mut cubic, MSS, Instant::from_millis(time));
        assert!(cubic.window() >= initial_cwnd);
    }

    #[test]
    fn fast_convergence_reduces_w_max() {
        let mut cubic = Cubic::new();
        cubic.set_mss(MSS);
        cubic.w_max = MSS * 50;
        cubic.cwnd = MSS * 30;

        // Loss while cwnd < w_max (a new competing flow) should pull w_max down.
        let w_max_prev = cubic.w_max;
        cubic.on_loss(Instant::from_millis(0), cubic.cwnd);
        assert!(cubic.w_max < w_max_prev);
    }

    #[test]
    fn test_cube_root() {
        for n in (1..1000000).step_by(99) {
            let a = n as f64;
            let a = a * a * a;
            let result = cube_root(a);
            println!("cube_root({a}) = {}", result.unwrap());
        }
    }

    #[test]
    #[should_panic]
    fn cube_root_zero() {
        cube_root(0.0).unwrap();
    }
}
