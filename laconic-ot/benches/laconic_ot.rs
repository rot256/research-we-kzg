use ark_bls12_381::{Bls12_381, Fr};
use ark_poly::Radix2EvaluationDomain;
use ark_std::rand::Rng;
use ark_std::test_rng;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use laconic_ot::{Choice, CommitmentKey, LaconicOTRecv, LaconicOTSender};

const MIN_LOG_SIZE: usize = 5;
const MAX_LOG_SIZE: usize = 19;

fn laconic_ot_benchmarks(c: &mut Criterion) {
    let name = "laconic_ot";
    let mut commit_benchmarks = c.benchmark_group(format!("{0}/commit", name));

    commit_benchmarks.sample_size(10);

    for log_len in MIN_LOG_SIZE..MAX_LOG_SIZE {
        commit_benchmarks.bench_with_input(
            BenchmarkId::from_parameter(log_len),
            &log_len,
            |b, _| {
                let rng = &mut test_rng();
                let num = 1 << log_len;

                let mut bits = Vec::with_capacity(log_len);
                for _ in 0..num {
                    bits.push(Choice::random(rng));
                }

                b.iter(|| {
                    let ck =
                        CommitmentKey::<Bls12_381, Radix2EvaluationDomain<Fr>>::setup(rng, num)
                            .unwrap();

                    let _sender = LaconicOTRecv::new(&ck, &bits);
                })
            },
        );
    }
    commit_benchmarks.finish();

    let mut send_benchmarks = c.benchmark_group(format!("{0}/send", name));
    for log_len in MIN_LOG_SIZE..MAX_LOG_SIZE {
        send_benchmarks.bench_with_input(BenchmarkId::from_parameter(log_len), &log_len, |b, _| {
            let rng = &mut test_rng();
            let num = 1 << log_len;

            let mut bits = Vec::with_capacity(log_len);
            for _ in 0..num {
                bits.push(Choice::random(rng));
            }

            let ck =
                CommitmentKey::<Bls12_381, Radix2EvaluationDomain<Fr>>::setup(rng, num).unwrap();
            let recv = LaconicOTRecv::new(&ck, &bits);

            let m0 = [0u8; 32];
            let m1 = [1u8; 32];

            b.iter(|| {
                let i = rng.gen_range(0..num);
                let sender = LaconicOTSender::new(&ck, recv.commitment());
                let _msg = sender.send(rng, i, m0, m1);
            })
        });
    }
    send_benchmarks.finish();

    let mut recv_benchmarks = c.benchmark_group(format!("{0}/recv", name));

    for log_len in MIN_LOG_SIZE..MAX_LOG_SIZE {
        recv_benchmarks.bench_with_input(BenchmarkId::from_parameter(log_len), &log_len, |b, _| {
            let rng = &mut test_rng();
            let num = 1 << log_len;

            let mut bits = Vec::with_capacity(log_len);
            for _ in 0..num {
                bits.push(Choice::random(rng));
            }

            let ck =
                CommitmentKey::<Bls12_381, Radix2EvaluationDomain<Fr>>::setup(rng, num).unwrap();
            let recv = LaconicOTRecv::new(&ck, &bits);

            let m0 = [0u8; 32];
            let m1 = [1u8; 32];

            let sender = LaconicOTSender::new(&ck, recv.commitment());

            let i = rng.gen_range(0..num);
            let msg = sender.send(rng, i, m0, m1);

            b.iter(|| {
                let _res = recv.recv(i, msg.clone());
            })
        });
    }
}

criterion_group! {
    name = laconic_ot;
    config = Criterion::default().sample_size(10); // .nresamples(10);
    targets = laconic_ot_benchmarks, // ipa_benchmarks
}
criterion_main!(laconic_ot);
