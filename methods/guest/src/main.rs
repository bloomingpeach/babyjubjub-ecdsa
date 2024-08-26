#![no_main]
#![no_std]

use alloc::string::String;
use alloc::vec::Vec;
use core::convert::TryInto;
use hyle_contract::{HyleInput, HyleOutput};
use hyle_guest::babyjubjub::{decompress_point, decompress_signature, verify};
use num_bigint::BigInt;
use risc0_zkvm::guest::env;
use serde::{Deserialize, Deserializer};

extern crate alloc;

risc0_zkvm::guest::entry!(ecdsa_verify);

#[derive(Deserialize)]
struct ECDSAVerifierInput {
    public_key: [u8; 32],
    message: String,
    #[serde(deserialize_with = "deserialize_fixed_size_array")]
    compressed_signature: [u8; 64],
}

fn deserialize_fixed_size_array<'de, D>(deserializer: D) -> Result<[u8; 64], D::Error>
where
    D: Deserializer<'de>,
{
    let vec: Vec<u8> = Vec::deserialize(deserializer)?;
    vec.try_into()
        .map_err(|_| serde::de::Error::custom("Expected 64 bytes"))
}

impl ECDSAVerifierInput {
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.public_key);
        bytes.extend_from_slice(self.message.as_bytes());
        bytes.extend_from_slice(&self.compressed_signature);
        bytes
    }
}

pub fn ecdsa_verify() {
    let input: HyleInput<ECDSAVerifierInput> = env::read();
    let public_key = input.program_inputs.public_key;
    let message = input.program_inputs.message;
    let signature = input.program_inputs.compressed_signature;

    let decompressed_public_key = decompress_point(public_key).expect("Invalid public key");
    let decompressed_signature = decompress_signature(&signature).expect("Invalid signature");

    // Convert message to BigInt
    let msg = BigInt::from_bytes_be(num_bigint::Sign::Plus, message.as_bytes());

    let is_valid = verify(decompressed_public_key, decompressed_signature, msg);

    env::commit(&HyleOutput {
        version: 1,
        identity: input.identity,
        tx_hash: input.tx_hash,
        program_outputs: if is_valid {
            "Valid signature"
        } else {
            "Invalid signature"
        },
        initial_state: input.initial_state,
        next_state: if is_valid {
            alloc::vec![1]
        } else {
            alloc::vec![0]
        },
    })
}
