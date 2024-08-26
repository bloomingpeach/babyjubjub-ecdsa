// BabyJubJub elliptic curve implementation in Rust.
// For LICENSE check https://github.com/arnaucube/babyjubjub-rs

use poseidon_rs::Poseidon;
pub type Fr = poseidon_rs::Fr; // alias
use crate::utils;
use arrayref::array_ref;

#[cfg(not(feature = "aarch64"))]
use blake_hash::Digest; // compatible version with Blake used at circomlib

#[cfg(feature = "aarch64")]
extern crate blake; // compatible version with Blake used at circomlib

use std::cmp::min;

use num_bigint::{BigInt, Sign, ToBigInt};
use num_traits::One;

use generic_array::GenericArray;

use lazy_static::lazy_static;

lazy_static! {
    static ref D: Fr = Fr::from_str("168696").unwrap();
    static ref D_BIG: BigInt = BigInt::parse_bytes(b"168696", 10).unwrap();
    static ref A: Fr = Fr::from_str("168700").unwrap();
    static ref A_BIG: BigInt = BigInt::parse_bytes(b"168700", 10).unwrap();
    pub static ref Q: BigInt = BigInt::parse_bytes(
        b"21888242871839275222246405745257275088548364400416034343698204186575808495617",10
    )
        .unwrap();
    static ref B8: Point = Point {
        x: Fr::from_str(
               "5299619240641551281634865583518297030282874472190772894086521144482721001553",
           )
            .unwrap(),
            y: Fr::from_str(
                "16950150798460657717958625567821834550301663161624707787222815936182638968203",
            )
                .unwrap(),
    };
    static ref ORDER: Fr = Fr::from_str(
        "21888242871839275222246405745257275088614511777268538073601725287587578984328",
    )
        .unwrap();

    // SUBORDER = ORDER >> 3
    static ref SUBORDER: BigInt = &BigInt::parse_bytes(
        b"21888242871839275222246405745257275088614511777268538073601725287587578984328",
        10,
    )
        .unwrap()
        >> 3;
    static ref POSEIDON: poseidon_rs::Poseidon = Poseidon::new();
}

#[derive(Clone, Debug)]
pub struct PointProjective {
    pub x: Fr,
    pub y: Fr,
    pub z: Fr,
}

impl PointProjective {
    pub fn affine(&self) -> Point {
        if self.z.is_zero() {
            return Point {
                x: Fr::zero(),
                y: Fr::zero(),
            };
        }

        let zinv = self.z.inverse().unwrap();
        let mut x = self.x;
        x.mul_assign(&zinv);
        let mut y = self.y;
        y.mul_assign(&zinv);

        Point { x, y }
    }

    #[allow(clippy::many_single_char_names)]
    pub fn add(&self, q: &PointProjective) -> PointProjective {
        // add-2008-bbjlp https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#addition-add-2008-bbjlp
        let mut a = self.z;
        a.mul_assign(&q.z);
        let mut b = a;
        b.square();
        let mut c = self.x;
        c.mul_assign(&q.x);
        let mut d = self.y;
        d.mul_assign(&q.y);
        let mut e = *D;
        e.mul_assign(&c);
        e.mul_assign(&d);
        let mut f = b;
        f.sub_assign(&e);
        let mut g = b;
        g.add_assign(&e);
        let mut x1y1 = self.x;
        x1y1.add_assign(&self.y);
        let mut x2y2 = q.x;
        x2y2.add_assign(&q.y);
        let mut aux = x1y1;
        aux.mul_assign(&x2y2);
        aux.sub_assign(&c);
        aux.sub_assign(&d);
        let mut x3 = a;
        x3.mul_assign(&f);
        x3.mul_assign(&aux);
        let mut ac = *A;
        ac.mul_assign(&c);
        let mut dac = d;
        dac.sub_assign(&ac);
        let mut y3 = a;
        y3.mul_assign(&g);
        y3.mul_assign(&dac);
        let mut z3 = f;
        z3.mul_assign(&g);

        PointProjective {
            x: x3,
            y: y3,
            z: z3,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Point {
    pub x: Fr,
    pub y: Fr,
}

impl Point {
    pub fn projective(&self) -> PointProjective {
        PointProjective {
            x: self.x,
            y: self.y,
            z: Fr::one(),
        }
    }

    pub fn mul_scalar(&self, n: &BigInt) -> Point {
        let mut r: PointProjective = PointProjective {
            x: Fr::zero(),
            y: Fr::one(),
            z: Fr::one(),
        };
        let mut exp: PointProjective = self.projective();
        let (_, b) = n.to_bytes_le();
        for i in 0..n.bits() {
            if test_bit(&b, i.try_into().unwrap()) {
                r = r.add(&exp);
            }
            exp = exp.add(&exp);
        }
        r.affine()
    }

    pub fn compress(&self) -> [u8; 32] {
        let p = &self;
        let mut r: [u8; 32] = [0; 32];
        let x_big = BigInt::parse_bytes(to_hex(&p.x).as_bytes(), 16).unwrap();
        let y_big = BigInt::parse_bytes(to_hex(&p.y).as_bytes(), 16).unwrap();
        let (_, y_bytes) = y_big.to_bytes_le();
        let len = min(y_bytes.len(), r.len());
        r[..len].copy_from_slice(&y_bytes[..len]);
        if x_big > (&Q.clone() >> 1) {
            r[31] |= 0x80;
        }
        r
    }

    pub fn equals(&self, p: Point) -> bool {
        if self.x == p.x && self.y == p.y {
            return true;
        }
        false
    }
}

pub fn test_bit(b: &[u8], i: usize) -> bool {
    b[i / 8] & (1 << (i % 8)) != 0
}

pub fn decompress_point(bb: [u8; 32]) -> Result<Point, String> {
    // https://tools.ietf.org/html/rfc8032#section-5.2.3
    let mut sign: bool = false;
    let mut b = bb;
    if b[31] & 0x80 != 0x00 {
        sign = true;
        b[31] &= 0x7F;
    }
    let y: BigInt = BigInt::from_bytes_le(Sign::Plus, &b[..]);
    if y >= Q.clone() {
        return Err("y outside the Finite Field over R".to_string());
    }
    let one: BigInt = One::one();

    // x^2 = (1 - y^2) / (a - d * y^2) (mod p)
    let den = utils::modinv(
        &utils::modulus(
            &(&A_BIG.clone() - utils::modulus(&(&D_BIG.clone() * (&y * &y)), &Q)),
            &Q,
        ),
        &Q,
    )?;
    let mut x: BigInt = utils::modulus(&((one - utils::modulus(&(&y * &y), &Q)) * den), &Q);
    x = utils::modsqrt(&x, &Q)?;

    if sign && (x <= (&Q.clone() >> 1)) || (!sign && (x > (&Q.clone() >> 1))) {
        x *= -(1.to_bigint().unwrap());
    }
    x = utils::modulus(&x, &Q);
    let x_fr: Fr = Fr::from_str(&x.to_string()).unwrap();
    let y_fr: Fr = Fr::from_str(&y.to_string()).unwrap();
    Ok(Point { x: x_fr, y: y_fr })
}

#[cfg(not(feature = "aarch64"))]
fn blh(b: &[u8]) -> Vec<u8> {
    let hash = blake_hash::Blake512::digest(b);
    hash.to_vec()
}

#[cfg(feature = "aarch64")]
fn blh(b: &[u8]) -> Vec<u8> {
    let mut hash = [0; 64];
    blake::hash(512, b, &mut hash).unwrap();
    hash.to_vec()
}

#[derive(Debug, Clone)]
pub struct Signature {
    pub r_b8: Point,
    pub s: BigInt,
}

impl Signature {
    pub fn compress(&self) -> [u8; 64] {
        let mut b: Vec<u8> = Vec::new();
        b.append(&mut self.r_b8.compress().to_vec());
        let (_, s_bytes) = self.s.to_bytes_le();
        let mut s_32bytes: [u8; 32] = [0; 32];
        let len = min(s_bytes.len(), s_32bytes.len());
        s_32bytes[..len].copy_from_slice(&s_bytes[..len]);
        b.append(&mut s_32bytes.to_vec());
        let mut r: [u8; 64] = [0; 64];
        r[..].copy_from_slice(&b[..]);
        r
    }
}

pub fn decompress_signature(b: &[u8; 64]) -> Result<Signature, String> {
    let r_b8_bytes: [u8; 32] = *array_ref!(b[..32], 0, 32);
    let s: BigInt = BigInt::from_bytes_le(Sign::Plus, &b[32..]);
    let r_b8 = decompress_point(r_b8_bytes);
    match r_b8 {
        Result::Err(err) => Err(err),
        Result::Ok(res) => Ok(Signature { r_b8: res, s }),
    }
}

pub struct PrivateKey {
    pub key: [u8; 32],
}

impl PrivateKey {
    pub fn import(b: Vec<u8>) -> Result<PrivateKey, String> {
        if b.len() != 32 {
            return Err(String::from("imported key can not be bigger than 32 bytes"));
        }
        let mut sk: [u8; 32] = [0; 32];
        sk.copy_from_slice(&b[..32]);
        Ok(PrivateKey { key: sk })
    }

    pub fn scalar_key(&self) -> BigInt {
        // not-compatible with circomlib implementation, but using Blake2b
        // let mut hasher = Blake2b::new();
        // hasher.update(sk_raw_bytes);
        // let mut h = hasher.finalize();

        // compatible with circomlib implementation
        let hash: Vec<u8> = blh(&self.key);
        let mut h: Vec<u8> = hash[..32].to_vec();

        // prune buffer following RFC 8032
        // https://tools.ietf.org/html/rfc8032#page-13
        h[0] &= 0xF8;
        h[31] &= 0x7F;
        h[31] |= 0x40;

        let sk = BigInt::from_bytes_le(Sign::Plus, &h[..]);
        sk >> 3
    }

    pub fn public(&self) -> Point {
        B8.mul_scalar(&self.scalar_key())
    }

    pub fn sign(&self, msg: BigInt) -> Result<Signature, String> {
        if msg > Q.clone() {
            return Err("msg outside the Finite Field".to_string());
        }
        // let (_, sk_bytes) = self.key.to_bytes_le();
        // let mut hasher = Blake2b::new();
        // hasher.update(sk_bytes);
        // let mut h = hasher.finalize(); // h: hash(sk), s: h[32:64]
        let mut h: Vec<u8> = blh(&self.key);

        let (_, msg_bytes) = msg.to_bytes_le();
        let mut msg32: [u8; 32] = [0; 32];
        msg32[..msg_bytes.len()].copy_from_slice(&msg_bytes[..]);
        let msg_fr: Fr = Fr::from_str(&msg.to_string()).unwrap();

        // https://tools.ietf.org/html/rfc8032#section-5.1.6
        let s = GenericArray::<u8, generic_array::typenum::U32>::from_mut_slice(&mut h[32..64]);
        let r_bytes = utils::concatenate_arrays(s, &msg32);
        let r_hashed: Vec<u8> = blh(&r_bytes);
        let mut r = BigInt::from_bytes_le(Sign::Plus, &r_hashed[..]);
        r = utils::modulus(&r, &SUBORDER);
        let r_b8: Point = B8.mul_scalar(&r);
        let a = &self.public();

        let hm_input = vec![r_b8.x, r_b8.y, a.x, a.y, msg_fr];
        let hm = POSEIDON.hash(hm_input)?;

        let mut s = &self.scalar_key() << 3;
        let hm_b = BigInt::parse_bytes(to_hex(&hm).as_bytes(), 16).unwrap();
        s = hm_b * s;
        s = r + s;
        s %= &SUBORDER.clone();

        Ok(Signature { r_b8, s })
    }
}

pub fn schnorr_hash(pk: &Point, msg: BigInt, c: &Point) -> Result<BigInt, String> {
    if msg > Q.clone() {
        return Err("msg outside the Finite Field".to_string());
    }
    let msg_fr: Fr = Fr::from_str(&msg.to_string()).unwrap();
    let hm_input = vec![pk.x, pk.y, c.x, c.y, msg_fr];
    let h = POSEIDON.hash(hm_input)?;
    let h_b = BigInt::parse_bytes(to_hex(&h).as_bytes(), 16).unwrap();
    Ok(h_b)
}

pub fn verify(pk: Point, sig: Signature, msg: BigInt) -> bool {
    if msg > Q.clone() {
        return false;
    }
    let msg_fr: Fr = Fr::from_str(&msg.to_string()).unwrap();
    let hm_input = vec![sig.r_b8.x, sig.r_b8.y, pk.x, pk.y, msg_fr];
    let hm = match POSEIDON.hash(hm_input) {
        Result::Err(_) => return false,
        Result::Ok(hm) => hm,
    };
    let l = B8.mul_scalar(&sig.s);
    let hm_b = BigInt::parse_bytes(to_hex(&hm).as_bytes(), 16).unwrap();
    let r = sig
        .r_b8
        .projective()
        .add(&pk.mul_scalar(&(8.to_bigint().unwrap() * hm_b)).projective());
    l.equals(r.affine())
}
