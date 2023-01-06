package ring

// NTT evaluates p2 = NTT(p1).
func (s *SubRing) NTT(p1, p2 []uint64) {
	s.Forward(s, p1, p2)
}

// NTTLazy evaluates p2 = NTT(p1) with p2 in [0, 2*modulus-1].
func (s *SubRing) NTTLazy(p1, p2 []uint64) {
	s.ForwardLazy(s, p1, p2)
}

// INTT evaluates p2 = INTT(p1).
func (s *SubRing) INTT(p1, p2 []uint64) {
	s.Backward(s, p1, p2)
}

// INTTLazy evaluates p2 = INTT(p1) with p2 in [0, 2*modulus-1].
func (s *SubRing) INTTLazy(p1, p2 []uint64) {
	s.BackwardLazy(s, p1, p2)
}
