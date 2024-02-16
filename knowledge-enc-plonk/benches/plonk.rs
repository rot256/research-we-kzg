// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! PLONK Benchmarks

use ark_bls12_377::Bls12_377;
use ark_bls12_381::Bls12_381;
use ark_ec::{PairingEngine, TEModelParameters};
use ark_ed_on_bls12_381::EdwardsParameters;
use ark_ff::{FftField, PrimeField};
use ark_std::test_rng;
use blake2;
use core::marker::PhantomData;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use plonk::commitment::{HomomorphicCommitment, IPA, KZG10};
use plonk::prelude::*;
use rand_core::OsRng;
use std::ops::Div;

/// Benchmark Circuit
#[derive(derivative::Derivative)]
#[derivative(Debug, Default)]
pub struct BenchCircuit<F, P> {
    /// Circuit Size
    size: usize,

    /// Field and parameters
    _phantom: PhantomData<(F, P)>,
}

impl<F, P> BenchCircuit<F, P> {
    /// Builds a new circuit with a constraint count of `2^degree`.
    #[inline]
    pub fn new(degree: usize) -> Self {
        Self {
            size: 1 << degree,
            _phantom: PhantomData::<(F, P)>,
        }
    }
}

impl<F, P> Circuit<F, P> for BenchCircuit<F, P>
where
    F: FftField + PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    const CIRCUIT_ID: [u8; 32] = [0xff; 32];

    #[inline]
    fn gadget(
        &mut self,
        composer: &mut StandardComposer<F, P>,
    ) -> Result<(), Error> {
        composer.add_dummy_lookup_table();
        while composer.circuit_bound() < self.size - 1 {
            composer.add_dummy_constraints();
        }
        Ok(())
    }

    #[inline]
    fn padded_circuit_size(&self) -> usize {
        self.size
    }
}

fn kzg10_benchmarks(c: &mut Criterion) {
    constraint_system_benchmark("KZG10", c);
}

/*
fn ipa_benchmarks(c: &mut Criterion) {
    constraint_system_benchmark::<
        <Bls12_377 as PairingEngine>::Fr,
        ark_ed_on_bls12_377::EdwardsParameters,
        IPA<<Bls12_377 as PairingEngine>::G1Affine, blake2::Blake2b>,
    >("IPA", c);
}
*/

/// Generates full benchmark suite for compiling, proving, and verifying.
fn constraint_system_benchmark(name: &str, c: &mut Criterion) {
    type E = Bls12_377;
    type P = ark_ed_on_bls12_377::EdwardsParameters;

    use ark_poly_commit::PolynomialCommitment;

    let label = b"ark".as_slice();

    const MINIMUM_DEGREE: usize = 5;
    const MAXIMUM_DEGREE: usize = 19;

    let pp = KZG10::<E>::setup(1 << MAXIMUM_DEGREE, None, &mut OsRng)
        .expect("Unable to sample public parameters.");

    let mut circuit = BenchCircuit::<<E as PairingEngine>::Fr, P>::new(10);
    let (pk_p, (vk, _pi_pos)) =
        circuit.compile(&pp).expect("Unable to compile circuit.");

    let (proof, pi, pck) = circuit.gen_proof(&pp, pk_p, &label).unwrap();

    let pvk =
        plonk::circuit::verify_proof::<<E as PairingEngine>::Fr, P, KZG10<E>>(
            &pp, vk, &proof, &pi, &label,
        )
        .expect("Unable to verify benchmark circuit.");

    let pk = pvk.pk();
    let sk = pck.sk();

    println!("pk-size: {}", pk.size());
    println!("sk-size: {}", sk.size());

    /*
    let mut compiling_benchmarks =
        c.benchmark_group(format!("{0}/compile", name));
    for degree in MINIMUM_DEGREE..MAXIMUM_DEGREE {
        let mut circuit =
            BenchCircuit::<<E as PairingEngine>::Fr, P>::new(degree);
        compiling_benchmarks.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    circuit
                        .compile::<KZG10<E>>(&pp)
                        .expect("Unable to compile circuit.")
                })
            },
        );
    }
    compiling_benchmarks.finish();
    */

    let mut keygen_benchmarks = c.benchmark_group("keygen");
    for degree in MINIMUM_DEGREE..MAXIMUM_DEGREE {
        let mut circuit =
            BenchCircuit::<<E as PairingEngine>::Fr, P>::new(degree);
        let (pk_p, _) = circuit
            .compile::<KZG10<E>>(&pp)
            .expect("Unable to compile circuit.");
        keygen_benchmarks.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    circuit
                        .gen_proof::<KZG10<E>>(&pp, pk_p.clone(), &label)
                        .unwrap()
                })
            },
        );
    }
    keygen_benchmarks.finish();

    let mut encapsulate_benchmarks = c.benchmark_group("encapsulate");
    for degree in MINIMUM_DEGREE..MAXIMUM_DEGREE {
        // already received proof/public key
        let mut circuit =
            BenchCircuit::<<E as PairingEngine>::Fr, P>::new(degree);
        let (pk_p, (vk, _pi_pos)) =
            circuit.compile(&pp).expect("Unable to compile circuit.");

        let (proof, pi, _) = circuit.gen_proof(&pp, pk_p, &label).unwrap();

        // verify the proof, producing a KEM key and encapsulate
        encapsulate_benchmarks.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    // verify the proof, producing a KEM key
                    let pvk = plonk::circuit::verify_proof::<
                        <E as PairingEngine>::Fr,
                        P,
                        KZG10<E>,
                    >(
                        &pp, vk.clone(), &proof, &pi, &label
                    )
                    .expect("Unable to verify benchmark circuit.");

                    // encapsulate to the proof
                    let pk = pvk.pk();
                    let rng = &mut test_rng();
                    let _ct = pk.enc(&pvk, rng);
                })
            },
        );
    }
    encapsulate_benchmarks.finish();

    let mut decapsulate_benchmarks = c.benchmark_group("decapsulate");
    for degree in MINIMUM_DEGREE..MAXIMUM_DEGREE {
        let mut circuit =
            BenchCircuit::<<E as PairingEngine>::Fr, P>::new(degree);
        let (pk_p, (vk, _pi_pos)) =
            circuit.compile(&pp).expect("Unable to compile circuit.");

        let (proof, pi, pck) = circuit.gen_proof(&pp, pk_p, &label).unwrap();

        let pvk = plonk::circuit::verify_proof::<
            <E as PairingEngine>::Fr,
            P,
            KZG10<E>,
        >(&pp, vk, &proof, &pi, &label)
        .expect("Unable to verify benchmark circuit.");

        let pk = pvk.pk();
        let sk = pck.sk();

        println!("pk-size: {}", pk.size());
        println!("sk-size: {}", sk.size());

        let rng = &mut test_rng();

        // encapsulate to the proof
        let (key1, ct) = pk.enc(&pvk, rng);

        decapsulate_benchmarks.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    let key2 = sk.dec(&pck, &ct);
                })
            },
        );
    }
    decapsulate_benchmarks.finish();
}

criterion_group! {
    name = plonk;
    config = Criterion::default().sample_size(10);
    targets = kzg10_benchmarks, // ipa_benchmarks
}
criterion_main!(plonk);
