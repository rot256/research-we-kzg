use std::{marker::PhantomData, ops::Div};

use crate::{CommitterKey, KZGKem, VerifierKey};
use ark_ec::AffineCurve;
use ark_ec::PairingEngine;
use ark_ff::ToBytes;
use ark_poly::UVPolynomial;
use ark_poly_commit::LabeledCommitment;
use ark_poly_commit::{kzg10, sonic_pc::SonicKZG10};
use ark_poly_commit::{LabeledPolynomial, PolynomialCommitment};
use ark_serialize::Write;
use ark_std::rand::Rng;
use ark_std::UniformRand;

type Key = Vec<u8>;

use blake2::Blake2s256;
use blake2::Digest;

fn to_key<E: PairingEngine>(m: E::Fqk) -> Key {
    let mut hsh = Blake2s256::new();
    m.write(&mut hsh).unwrap();
    hsh.finalize().to_vec()
}

pub struct Ct<E: PairingEngine, P: UVPolynomial<E::Fr, Point = E::Fr>> {
    _ph: PhantomData<(E, P)>,
    h: E::G2Affine,
}

pub struct Cts<E: PairingEngine, P: UVPolynomial<E::Fr, Point = E::Fr>> {
    _ph: PhantomData<(E, P)>,
    cts: Vec<Ct<E, P>>,
}

pub struct Claim<E: PairingEngine, P: UVPolynomial<E::Fr, Point = E::Fr>> {
    comm: kzg10::Commitment<E>,
    rand: kzg10::Randomness<E::Fr, P>,
    poly: P,
    pont: P::Point,
}

impl<E: PairingEngine, P: UVPolynomial<E::Fr, Point = E::Fr>> Claim<E, P> {
    fn open(
        &self,
        crs: &<KZGKem<E, P> as PolynomialCommitment<E::Fr, P>>::CommitterKey,
    ) -> KemSK<E, P>
    where
        for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    {
        use ark_std::One;
        let crs = &crs.inner;
        let poly = LabeledPolynomial::new("test".to_string(), self.poly.clone(), None, None);
        let comm = LabeledCommitment::new("test".to_string(), self.comm.clone(), None);
        let proof = SonicKZG10::open_individual_opening_challenges(
            &crs,
            vec![&poly],
            vec![&comm],
            &self.pont,
            &|_| E::Fr::one(),
            vec![&self.rand],
            None,
        )
        .unwrap();

        KemSK {
            open: proof.w,
            comm: self.comm.clone(),
            pont: self.pont,
            eval: self.poly.evaluate(&self.pont),
        }
    }
}

pub struct KemSK<E: PairingEngine, P: UVPolynomial<E::Fr, Point = E::Fr>> {
    open: E::G1Affine,
    comm: kzg10::Commitment<E>,
    pont: P::Point,
    eval: E::Fr,
}

impl<E: PairingEngine, P: UVPolynomial<E::Fr, Point = E::Fr>> KemSK<E, P> {
    pub fn random<R: Rng>(
        rng: &mut R,
        crs: &<KZGKem<E, P> as PolynomialCommitment<E::Fr, P>>::CommitterKey,
        deg: usize,
    ) -> Self
    where
        for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    {
        let poly = P::rand(deg, rng);
        let poly = LabeledPolynomial::new("test".to_string(), poly, None, None);
        let pont = E::Fr::rand(rng);
        let (comms, rands) = KZGKem::commit(crs, vec![&poly], Some(rng)).unwrap();
        let pont = poly.evaluate(&pont);
        Claim {
            comm: comms[0].commitment().clone(),
            rand: rands[0].clone(),
            poly: poly.polynomial().clone(),
            pont,
        }
        .open(crs)
    }

    pub fn pk(&self) -> KemPK<E, P> {
        KemPK {
            comm: self.comm.clone(),
            pont: self.pont,
            eval: self.eval,
        }
    }

    pub fn dec(
        &self,
        crs: &<KZGKem<E, P> as PolynomialCommitment<E::Fr, P>>::CommitterKey,
        ct: &Ct<E, P>,
    ) -> Key
    where
        for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    {
        let m = E::pairing(self.open, ct.h);
        to_key::<E>(m)
    }
}

pub struct KemPK<E: PairingEngine, P: UVPolynomial<E::Fr, Point = E::Fr>> {
    comm: kzg10::Commitment<E>,
    pont: P::Point,
    eval: E::Fr,
}

impl<E: PairingEngine, P: UVPolynomial<E::Fr, Point = E::Fr>> KemPK<E, P> {
    pub fn enc<R: Rng>(
        &self,
        crs: &<KZGKem<E, P> as PolynomialCommitment<E::Fr, P>>::VerifierKey,
        rng: &mut R,
    ) -> (Key, Ct<E, P>)
    where
        for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    {
        let crs = &crs.inner;

        let r = E::Fr::rand(rng);
        let x = self.pont;
        let y = self.eval;
        let c = self.comm.0;

        let l = c.mul(r) - crs.g.mul(r * y); // e * (c - [y])
        let h = crs.beta_h.mul(r) - crs.h.mul(x * r); // h = r * ([t] - [x])

        // key = e(l, [1])
        let m = E::pairing(l, crs.h);

        (
            to_key::<E>(m),
            Ct {
                _ph: PhantomData,
                h: h.into(),
            },
        )
    }
}

pub struct MultiSK<E: PairingEngine, P: UVPolynomial<E::Fr, Point = E::Fr>> {
    sks: Vec<KemSK<E, P>>,
}

impl<E: PairingEngine, P: UVPolynomial<E::Fr, Point = E::Fr>> MultiSK<E, P> {
    pub fn dec(
        &self,
        crs: &<KZGKem<E, P> as PolynomialCommitment<E::Fr, P>>::CommitterKey,
        cts: &Cts<E, P>,
    ) -> Key
    where
        for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    {
        let mut hsh = Blake2s256::new();
        for (sk, ct) in self.sks.iter().zip(cts.cts.iter()) {
            let key = sk.dec(crs, ct);
            hsh.write(&key).unwrap();
        }
        hsh.finalize().to_vec()
    }
}

impl<E: PairingEngine, P: UVPolynomial<E::Fr, Point = E::Fr>> MultiPK<E, P> {
    pub fn enc<R: Rng>(
        &self,
        crs: &<KZGKem<E, P> as PolynomialCommitment<E::Fr, P>>::VerifierKey,
        rng: &mut R,
    ) -> (Key, Cts<E, P>)
    where
        for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    {
        let mut hsh = Blake2s256::new();
        let mut cts = Vec::new();
        for pk in self.pks.iter() {
            let (key, ct) = pk.enc(crs, rng);
            hsh.write(&key).unwrap();
            cts.push(ct);
        }
        (
            hsh.finalize().to_vec(),
            Cts {
                _ph: PhantomData,
                cts,
            },
        )
    }
}

pub struct MultiPK<E: PairingEngine, P: UVPolynomial<E::Fr, Point = E::Fr>> {
    pks: Vec<KemPK<E, P>>,
}

impl<E, P> VerifierKey<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    pub fn pk(&self) -> MultiPK<E, P> {
        let mut pks = Vec::new();
        for (comm, x, y) in self.trace.lock().unwrap().elem.iter() {
            pks.push(KemPK {
                comm: comm.clone(),
                pont: x.clone(),
                eval: y.clone(),
            });
        }
        MultiPK { pks }
    }
}

impl<E, P> CommitterKey<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    pub fn sk(&self) -> MultiSK<E, P> {
        let mut sks = Vec::new();
        for (comm, rand, poly, pont) in self.trace.lock().unwrap().elem.iter() {
            sks.push(
                Claim {
                    comm: comm.clone(),
                    rand: rand.clone(),
                    poly: poly.clone(),
                    pont: pont.clone(),
                }
                .open(self),
            );
        }
        MultiSK { sks }
    }
}

#[test]
fn test_kem() {
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_poly::univariate::DensePolynomial;
    use ark_std::test_rng;

    let rng = &mut test_rng();
    let crs = KZGKem::<Bls12_381, DensePolynomial<Fr>>::setup(100, None, rng).unwrap();
    let (ck, vk) = KZGKem::<Bls12_381, DensePolynomial<Fr>>::trim(&crs, 100, 0, None).unwrap();

    let sk = KemSK::random(rng, &ck, 10);
    let pk = sk.pk();

    let (key1, ct) = pk.enc(&vk, rng);
    let key2 = sk.dec(&ck, &ct);

    assert_eq!(key1, key2);
}
