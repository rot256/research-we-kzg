use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_poly::EvaluationDomain;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::UniformRand;
use ark_std::Zero;
use std::marker::PhantomData;
use std::ops::Mul;

use crate::kzg_fk_open::precompute_y;

/// Simulation-Extractable vector commitment based on KZG
pub struct VcKZG<E: Pairing, D: EvaluationDomain<E::ScalarField>> {
    _e: PhantomData<E>,
    _d: PhantomData<D>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, PartialEq, Eq, Debug)]
pub struct CommitmentKey<E: Pairing, D: EvaluationDomain<E::ScalarField>> {
    /// length of messages to which we commit,
    /// This is called ell in the paper
    pub message_length: usize,

    /// evaluation domain that we use to represent vectors
    /// should support polynomials of degree deg >= ell+1
    pub domain: D,

    /// powers-of-alpha: u[i] = g1^{alpha^i}
    /// i should range from 0 to deg
    /// Note: u[0] = g1
    pub u: Vec<E::G1Affine>,

    // /// lagrange version of u
    // /// u_lag[i] = g1^{l_i(alpha)}, where
    // /// l_i is the ith lagrange polynomial
    // /// i should range from 0 to deg
    // pub u_lag: Vec<E::G1Affine>,
    /// same as u, but for the hiding part
    /// hat_u[i] = h1^{alpha^i}
    /// i should range from 0 to deg
    pub hat_u: Vec<E::G1Affine>,

    /// lagrange version of u and hat_u
    /// Let l_i be the ith lagrange poly. Then:
    /// lag[i] = g1^{l_i(alpha)}
    /// lag[deg+i] = h1^{l_i(alpha)}
    /// for i in 0..deg
    pub lagranges: Vec<E::G1Affine>,

    /// generator of G2
    pub g2: E::G2Affine,

    /// r = g2^{\alpha}, needed for verification
    pub r: E::G2Affine,

    /// precomputed denominators in the exponent
    /// all prepared for pairing
    /// namely, d[i] = g2^{alpha - zi},
    /// where zi is the ith evaluation point
    /// i should range from 0 to deg
    pub d: Vec<E::G2Affine>,

    /// y = DFT_{2d}(hat_s) for
    /// hat_s = [u[d-1],...,u[0], d+2 neutral elements]
    /// precomputed for use in the FK technique
    pub y: Vec<E::G1Affine>,
}

#[derive(CanonicalSerialize)]
pub struct Opening<E: Pairing> {
    /// commitment to witness polynomial g1^{psi(alpha)}
    pub v: E::G1Affine,
}

#[derive(CanonicalSerialize)]
pub struct Commitment<E: Pairing> {
    pub com_kzg: E::G1Affine,
}

pub struct State<E: Pairing> {
    /// stores both the evaluations of the polynomial
    /// and the evaluations of the masking polynomial
    /// polynomial: 0..deg, masking: deg..2*deg
    pub evals: Vec<E::ScalarField>,

    /// optionally stores precomputed KZG openings
    /// Note: this is only the group element part
    pub precomputed_v: Option<Vec<E::G1>>,
}

impl<E: Pairing, D: EvaluationDomain<E::ScalarField>> CommitmentKey<E, D> {
    pub fn setup<R: rand::Rng>(
        rng: &mut R,
        message_length: usize,
    ) -> Result<CommitmentKey<E, D>, ()> {
        if message_length < 1 {
            return Err(());
        }

        // generate an evaluation domain
        // should support polynomials to degree >= message_length + 1
        let domain = D::new(message_length);
        if domain.is_none() {
            return Err(());
        }
        let domain = domain.unwrap();

        // sample generators g1 and g2
        let g1 = E::G1::rand(rng);
        let g2 = E::G2::rand(rng);
        if g1.is_zero() || g2.is_zero() {
            return Err(());
        }

        // sample hiding generator h
        let h = E::G1::rand(rng);

        // sample secret exponent alpha
        let alpha = E::ScalarField::rand(rng);

        // raise g1 to the powers of alpha --> u
        // raise h to the powers of alpha  --> hat_u
        let deg = domain.size() - 1;
        let mut u: Vec<E::G1Affine> = Vec::new();
        let mut hat_u: Vec<E::G1Affine> = Vec::new();
        let mut curr_g = g1;
        let mut curr_h = h;
        u.push(curr_g.into_affine());
        hat_u.push(curr_h.into_affine());
        for _ in 1..=deg {
            curr_g = curr_g.mul(alpha);
            u.push(curr_g.into_affine());
            curr_h = curr_h.mul(alpha);
            hat_u.push(curr_h.into_affine());
        }

        // compute exponentiated lagrange coefficients
        // Note: If a standard powers-of-tau setup is used,
        // this can be publicly computed from u and hat_u
        let lf = domain.evaluate_all_lagrange_coefficients(alpha);
        let mut lagranges = Vec::with_capacity(2 * deg);
        for i in 0..=deg {
            lagranges.push(u[0].mul(lf[i]).into_affine());
        }

        //compute r = g2^{alpha}
        let r = g2.mul(alpha).into_affine();

        // compute all d[i] = g2^{alpha - zi}
        let mut d = Vec::new();
        for i in 0..message_length {
            let z = domain.element(i);
            let exponent: E::ScalarField = alpha - z;
            d.push(g2.mul(exponent).into_affine());
        }

        // precompute y and hat_y for FK algorithm
        let y = precompute_y::<E, D>(&u, &domain);

        // assemble commitment key
        let g2 = g2.into_affine();
        Ok(CommitmentKey {
            message_length,
            domain,
            u,
            hat_u,
            lagranges,
            g2,
            r,
            d,
            y,
        })
    }
}
