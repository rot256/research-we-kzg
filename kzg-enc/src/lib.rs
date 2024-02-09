use std::rc::Rc;
use std::sync::Mutex;
use std::{marker::PhantomData, ops::Div};

use ark_ec::AffineCurve;
use ark_ec::PairingEngine;
use ark_ff::ToBytes;
use ark_poly::{polynomial, Polynomial, UVPolynomial};
use ark_poly_commit::PCRandomness;
use ark_poly_commit::{kzg10, sonic_pc::SonicKZG10, PCProof, PolynomialCommitment};
use ark_poly_commit::{PCCommitterKey, PCPreparedVerifierKey, PCVerifierKey};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError, Write};
use ark_std::Zero;

mod kem;

/// This implements polynomial commit
///
/// When opening, it returns nothing.
/// When verifying it stores the commitment and the opening.
pub struct KZGKem<E: PairingEngine, P: Polynomial<E::Fr>> {
    _ph: PhantomData<(E, P)>,
}

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct NonProof;

impl PCProof for NonProof {
    fn size_in_bytes(&self) -> usize {
        0
    }
}

impl ToBytes for NonProof {
    fn write<W: ark_std::io::Write>(&self, _writer: W) -> Result<(), ark_std::io::Error> {
        Ok(())
    }
}

impl CanonicalDeserialize for NonProof {
    fn deserialize<R: ark_std::io::Read>(_reader: R) -> Result<Self, SerializationError> {
        Ok(Self)
    }
}

impl CanonicalSerialize for NonProof {
    fn serialized_size(&self) -> usize {
        0
    }

    fn serialize<W: Write>(&self, _writer: W) -> Result<(), SerializationError> {
        Ok(())
    }
}

impl Into<Vec<NonProof>> for NonProof {
    fn into(self) -> Vec<NonProof> {
        vec![self]
    }
}

impl From<Vec<NonProof>> for NonProof {
    fn from(_: Vec<NonProof>) -> Self {
        NonProof
    }
}

#[derive(Debug)]
struct TraceProver<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
{
    elem: Vec<(
        kzg10::Commitment<E>,
        kzg10::Randomness<E::Fr, P>,
        P,
        P::Point,
    )>,
}

#[derive(Debug)]
struct TraceVerifier<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
{
    elem: Vec<(kzg10::Commitment<E>, P::Point, E::Fr)>,
}

#[derive(Debug, Clone)]
pub struct CommitterKey<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    inner: <SonicKZG10<E, P> as PolynomialCommitment<E::Fr, P>>::CommitterKey,
    trace: Rc<Mutex<TraceProver<E, P>>>,
}

#[derive(Debug, Clone)]
pub struct VerifierKey<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    inner: <SonicKZG10<E, P> as PolynomialCommitment<E::Fr, P>>::VerifierKey,
    trace: Rc<Mutex<TraceVerifier<E, P>>>,
}

impl<E, P> PCVerifierKey for VerifierKey<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    fn max_degree(&self) -> usize {
        self.inner.max_degree()
    }

    fn supported_degree(&self) -> usize {
        self.inner.supported_degree()
    }
}

impl<E, P> CanonicalDeserialize for VerifierKey<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    fn deserialize<R: ark_std::io::Read>(_reader: R) -> Result<Self, SerializationError> {
        unimplemented!()
    }
}

impl<E, P> CanonicalSerialize for VerifierKey<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    fn serialized_size(&self) -> usize {
        unimplemented!()
    }

    fn serialize<W: Write>(&self, _writer: W) -> Result<(), SerializationError> {
        unimplemented!()
    }
}

impl<E, P> CanonicalDeserialize for CommitterKey<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    fn deserialize<R: ark_std::io::Read>(_reader: R) -> Result<Self, SerializationError> {
        unimplemented!()
    }
}

impl<E, P> CanonicalSerialize for CommitterKey<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    fn serialized_size(&self) -> usize {
        unimplemented!()
    }

    fn serialize<W: Write>(&self, _writer: W) -> Result<(), SerializationError> {
        unimplemented!()
    }
}

impl<E, P> PCCommitterKey for CommitterKey<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    fn max_degree(&self) -> usize {
        self.inner.max_degree()
    }

    fn supported_degree(&self) -> usize {
        self.inner.supported_degree()
    }
}

impl<E, P> PCPreparedVerifierKey<VerifierKey<E, P>> for VerifierKey<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    fn prepare(vk: &VerifierKey<E, P>) -> Self {
        vk.clone()
    }
}

impl<E, P> PolynomialCommitment<E::Fr, P> for KZGKem<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    type Commitment = kzg10::Commitment<E>;
    type UniversalParams = kzg10::UniversalParams<E>;
    type CommitterKey = CommitterKey<E, P>;
    type VerifierKey = VerifierKey<E, P>;
    type PreparedVerifierKey = VerifierKey<E, P>;
    type PreparedCommitment = kzg10::PreparedCommitment<E>;
    type Randomness = kzg10::Randomness<E::Fr, P>;

    // proofs are dummy
    type Proof = NonProof;
    type BatchProof = NonProof;

    type Error = <SonicKZG10<E, P> as PolynomialCommitment<E::Fr, P>>::Error;

    fn setup<R: ark_std::rand::prelude::RngCore>(
        max_degree: usize,
        num_vars: Option<usize>,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        SonicKZG10::<E, P>::setup(max_degree, num_vars, rng)
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        supported_hiding_bound: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        assert!(enforced_degree_bounds.is_none());

        let (ck, vk) = SonicKZG10::<E, P>::trim(
            pp,
            supported_degree,
            supported_hiding_bound,
            enforced_degree_bounds,
        )?;

        Ok((
            CommitterKey {
                inner: ck,
                trace: Rc::new(Mutex::new(TraceProver { elem: vec![] })),
            },
            VerifierKey {
                inner: vk,
                trace: Rc::new(Mutex::new(TraceVerifier { elem: vec![] })),
            },
        ))
    }

    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a ark_poly_commit::LabeledPolynomial<E::Fr, P>>,
        rng: Option<&mut dyn ark_std::rand::prelude::RngCore>,
    ) -> Result<
        (
            Vec<ark_poly_commit::LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    >
    where
        P: 'a,
    {
        SonicKZG10::<E, P>::commit(&ck.inner, polynomials, rng)
    }

    fn open_individual_opening_challenges<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a ark_poly_commit::LabeledPolynomial<E::Fr, P>>,
        commitments: impl IntoIterator<Item = &'a ark_poly_commit::LabeledCommitment<Self::Commitment>>,
        point: &'a <P as Polynomial<E::Fr>>::Point,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>, // randomness from kzg
        _rng: Option<&mut dyn ark_std::rand::prelude::RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a,
    {
        // compute random linear combination

        let mut poly_sum = P::zero();
        let mut comm_sum = E::G1Projective::zero();
        let mut rand_sum = kzg10::Randomness::empty();

        for (i, ((poly, comm), rand)) in labeled_polynomials
            .into_iter()
            .zip(commitments)
            .zip(rands.into_iter())
            .enumerate()
        {
            assert!(poly.degree_bound().is_none());
            assert!(poly.hiding_bound().is_none());

            let chal = opening_challenges(i as u64);
            comm_sum += comm.commitment().0.mul(chal);
            poly_sum += (chal, poly.polynomial());
            rand_sum += (chal, rand);
        }

        // store (comm_sum, rand_sum, poly_sum, point)
        ck.trace.lock().unwrap().elem.push((
            kzg10::Commitment(comm_sum.into()),
            rand_sum,
            poly_sum,
            *point,
        ));

        Ok(NonProof)
    }

    fn check_individual_opening_challenges<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a ark_poly_commit::LabeledCommitment<Self::Commitment>>,
        point: &'a <P as Polynomial<E::Fr>>::Point,
        values: impl IntoIterator<Item = E::Fr>,
        _proof: &Self::Proof, // the proof is dummy
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        _rng: Option<&mut dyn ark_std::rand::prelude::RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        // compute random linear combination

        let mut value_sum = E::Fr::zero();
        let mut comm_sum = E::G1Projective::zero();

        for (i, (value, comm)) in values.into_iter().zip(commitments.into_iter()).enumerate() {
            let chal = opening_challenges(i as u64);
            comm_sum += comm.commitment().0.mul(chal);
            value_sum += chal * value;
        }

        // store (comm_sum, point, values)
        vk.trace
            .lock()
            .unwrap()
            .elem
            .push((kzg10::Commitment(comm_sum.into()), *point, value_sum));

        Ok(true)
    }
}
