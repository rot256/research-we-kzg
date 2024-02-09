use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_poly::EvaluationDomain;
use ark_std::UniformRand;
use ark_std::Zero;

use std::ops::Mul;

use crate::kzg_types::Commitment;
use crate::kzg_types::Opening;
use crate::kzg_types::State;
use crate::kzg_utils::plain_kzg_com;
use crate::kzg_utils::witness_evals_inside;
use crate::{
    kzg_fk_open::precompute_y,
    kzg_types::{CommitmentKey, VcKZG},
};

impl<E: Pairing, D: EvaluationDomain<E::ScalarField>> VcKZG<E, D> {
    pub fn setup<R: rand::Rng>(
        rng: &mut R,
        message_length: usize,
    ) -> Result<CommitmentKey<E, D>, ()> {
        CommitmentKey::setup(rng, message_length)
    }

    pub fn commit<R: rand::Rng>(
        rng: &mut R,
        ck: &CommitmentKey<E, D>,
        m: &Vec<E::ScalarField>,
    ) -> (Commitment<E>, State<E>) {
        // evals[0..domain.size] will store evaluations of our polynomial
        // over our evaluation domain, namely
        // evals[i] = m[i]   if m[i] is defined,
        // evals[i] = random if not
        // evals[domain.size()..2*domain.size()] will store evaluations of
        // the random masking polynomial used for hiding
        // we keep both evaluations in the same vector so that
        // we can easily do a single MSM later
        let dsize = ck.domain.size();
        let mut evals = Vec::with_capacity(2 * dsize);
        for i in 0..m.len() {
            evals.push(m[i]);
        }
        for _ in m.len()..2 * ck.domain.size() {
            evals.push(E::ScalarField::rand(rng));
        }

        // from our evaluations, we compute a standard KZG commitment
        let com_kzg = plain_kzg_com(ck, &evals);

        let state = State {
            evals,
            precomputed_v: None,
        };
        let com = Commitment { com_kzg };
        (com, state)
    }

    pub fn open(ck: &CommitmentKey<E, D>, st: &State<E>, i: u32) -> Result<Opening<E>, ()> {
        if i as usize >= ck.message_length {
            return Err(());
        }

        // compute v: the KZG opening, which is a KZG commitment
        // to the witness polynomial. Either we already have it
        // precomputed, or we compute it in evaluation form
        let v = if let Some(vs) = &st.precomputed_v {
            vs[i as usize].into_affine()
        } else {
            let deg = ck.domain.size();
            let mut witn_evals = Vec::new();
            witness_evals_inside::<E, D>(&ck.domain, &st.evals, i as usize, &mut witn_evals);
            witness_evals_inside::<E, D>(
                &ck.domain,
                &st.evals[deg..2 * deg],
                i as usize,
                &mut witn_evals,
            );
            plain_kzg_com(&ck, &witn_evals)
        };

        Ok(Opening { v })
    }
}
