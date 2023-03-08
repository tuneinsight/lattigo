package advanced

import (
	"encoding/json"
)

// MarshalBinary returns a JSON reprsentation of the the target HomomorphicDFTMatrixLiteral on a slice of bytes.
// See `Marshal` from the `encoding/json` package.
func (d *HomomorphicDFTMatrixLiteral) MarshalBinary() (data []byte, err error) {
	return json.Marshal(d)
}

// UnmarshalBinary reads a JSON representation on the target HomomorphicDFTMatrixLiteral struct.
// See `Unmarshal` from the `encoding/json` package.
func (d *HomomorphicDFTMatrixLiteral) UnmarshalBinary(data []byte) error {
	return json.Unmarshal(data, d)
}

// MarshalBinary returns a JSON reprsentation of the the target EvalModLiteral struct on a slice of bytes.
// See `Marshal` from the `encoding/json` package.
func (evm *EvalModLiteral) MarshalBinary() (data []byte, err error) {
	return json.Marshal(evm)
}

// UnmarshalBinary reads a JSON representation on the target EvalModLiteral struct.
// See `Unmarshal` from the `encoding/json` package.
func (evm *EvalModLiteral) UnmarshalBinary(data []byte) (err error) {
	return json.Unmarshal(data, evm)
}
