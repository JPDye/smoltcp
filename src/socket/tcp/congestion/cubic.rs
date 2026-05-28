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

    #[test]
    fn test_cubic() {
        let remote_window = 64 * 1024 * 1024;
        let now = Instant::from_millis(0);

        for i in 0..10 {
            for j in 0..9 {
                let mut cubic = Cubic::new();
                // Set remote window.
                cubic.set_remote_window(remote_window);

                cubic.set_mss(1480);

                if i & 1 == 0 {
                    cubic.on_rto(now, cubic.window());
                } else {
                    cubic.on_dup_ack(now, 1480, cubic.window());
                }

                cubic.pre_transmit(now);

                let mut n = i;
                for _ in 0..j {
                    n *= i;
                }

                let elapsed = Instant::from_millis(n);
                cubic.pre_transmit(elapsed);

                let cwnd = cubic.window();
                println!("Cubic: elapsed = {}, cwnd = {}", elapsed, cwnd);

                assert!(cwnd >= cubic.min_cwnd);
                assert!(cubic.window() <= remote_window);
            }
        }
    }

    #[test]
    fn cubic_time_inversion() {
        let mut cubic = Cubic::new();

        let t1 = Instant::from_micros(0);
        let t2 = Instant::from_micros(i64::MAX);

        cubic.on_rto(t2, cubic.window());
        cubic.pre_transmit(t1);

        let cwnd = cubic.window();
        println!("Cubic:time_inversion: cwnd: {}, cubic: {cubic:?}", cwnd);

        assert!(cwnd >= cubic.min_cwnd);
        assert!(cwnd <= cubic.rwnd);
    }

    #[test]
    fn cubic_long_elapsed_time() {
        let mut cubic = Cubic::new();

        let t1 = Instant::from_millis(0);
        let t2 = Instant::from_micros(i64::MAX);

        cubic.on_rto(t1, cubic.window());
        cubic.pre_transmit(t2);

        let cwnd = cubic.window();
        println!("Cubic:long_elapsed_time: cwnd: {}", cwnd);

        assert!(cwnd >= cubic.min_cwnd);
        assert!(cwnd <= cubic.rwnd);
    }

    #[test]
    fn cubic_last_update() {
        let mut cubic = Cubic::new();

        let t1 = Instant::from_millis(0);
        let t2 = Instant::from_millis(100);
        let t3 = Instant::from_millis(199);
        let t4 = Instant::from_millis(20000);

        cubic.on_rto(t1, cubic.window());

        cubic.pre_transmit(t2);
        let cwnd2 = cubic.window();

        cubic.pre_transmit(t3);
        let cwnd3 = cubic.window();

        cubic.pre_transmit(t4);
        let cwnd4 = cubic.window();

        println!(
            "Cubic:last_update: cwnd2: {}, cwnd3: {}, cwnd4: {}",
            cwnd2, cwnd3, cwnd4
        );

        assert_eq!(cwnd2, cwnd3);
        assert_ne!(cwnd2, cwnd4);
    }

    #[test]
    fn cubic_slow_start() {
        let mut cubic = Cubic::new();

        let t1 = Instant::from_micros(0);

        let cwnd = cubic.window();
        let ack_len = 1024;

        cubic.on_ack(t1, ack_len, cubic.window(), &RttEstimator::default());

        assert!(cubic.window() > cwnd);

        for i in 1..1000 {
            let t2 = Instant::from_micros(i);
            cubic.on_ack(t2, ack_len * 100, cubic.window(), &RttEstimator::default());
            assert!(cubic.window() <= cubic.rwnd);
        }

        let t3 = Instant::from_micros(2000);

        let cwnd = cubic.window();
        cubic.on_rto(t3, cubic.window());
        assert_eq!(cwnd >> 1, cubic.ssthresh);
    }

    #[test]
    fn cubic_pre_transmit() {
        let mut cubic = Cubic::new();
        cubic.pre_transmit(Instant::from_micros(2000));
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
