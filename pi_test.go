package main

import (
	"fmt"
	"math/big"
	"testing"
)

// piRef is "3" followed by the first 1000+ decimal digits of π. Cross-checked
// against angio.net/pi and independent sources. Only the first 1000 decimals
// are relied upon (the trailing characters are kept for headroom but not used
// in assertions).
const piRef = "3." +
	"1415926535897932384626433832795028841971693993751058209749" +
	"4459230781640628620899862803482534211706798214808651328230" +
	"6647093844609550582231725359408128481117450284102701938521" +
	"1055596446229489549303819644288109756659334461284756482337" +
	"8678316527120190914564856692346034861045432664821339360726" +
	"0249141273724587006606315588174881520920962829254091715364" +
	"3678925903600113305305488204665213841469519415116094330572" +
	"7036575959195309218611738193261179310511854807446237996274" +
	"9567351885752724891227938183011949129833673362440656643086" +
	"0213949463952247371907021798609437027705392171762931767523" +
	"8467481846766940513200056812714526356082778577134275778960" +
	"9173637178721468440901224953430146549585371050792279689258" +
	"9235420199561121290219608640344181598136297747713099605187" +
	"0721134999999837297804995105973173281609631859502445945534" +
	"6908302642522308253344685035261931188171010003137838752886" +
	"5875332083814206171776691473035982534904287554687311595628" +
	"6388235378759375195778185778053217122680661300192787661119" +
	"5909216420198938095257201065485863278"

// refDigit returns the expected digit at digitPos using the program's 1-based
// convention: position 1 is the integer part '3', position N (N≥2) is the
// (N−1)th decimal digit.
func refDigit(digitPos int) int {
	if digitPos == 1 {
		return int(piRef[0] - '0') // '3'
	}
	return int(piRef[digitPos] - '0') // skip the '.' at index 1
}

func TestPiDecimalString(t *testing.T) {
	for _, d := range []int{50, 100, 250, 500, 1000} {
		got := piDecimalString(d)
		want := "3" + piRef[2:2+d]
		if got != want {
			// Report the first divergence to make transcription/algorithm bugs obvious.
			for i := 0; i < len(want) && i < len(got); i++ {
				if got[i] != want[i] {
					t.Fatalf("piDecimalString(%d) diverges at index %d: got %q want %q\n got=%s\nwant=%s",
						d, i, got[i], want[i], got, want)
				}
			}
			t.Fatalf("piDecimalString(%d): length %d, want %d", d, len(got), len(want))
		}
	}
}

func TestExtractDigit(t *testing.T) {
	positions := []int{1, 2, 3, 4, 5, 10, 50, 100, 250, 500, 762, 763, 764, 765, 766, 767, 768, 769, 999, 1000, 1001}
	for _, p := range positions {
		if got, want := extractDigit(p), refDigit(p); got != want {
			t.Errorf("extractDigit(%d) = %d, want %d", p, got, want)
		}
	}
}

// TestFeynmanPoint guards the six consecutive 9s at decimals 762–767 (positions
// 763–768) and the following 8 — the classic carry-propagation hazard.
func TestFeynmanPoint(t *testing.T) {
	for p := 763; p <= 768; p++ {
		if got := extractDigit(p); got != 9 {
			t.Errorf("extractDigit(%d) = %d, want 9 (Feynman point)", p, got)
		}
	}
	if got := extractDigit(769); got != 8 {
		t.Errorf("extractDigit(769) = %d, want 8", got)
	}
}

// TestRegressionLocks pins the documented examples so any optimization that
// changes a deep digit fails loudly.
func TestRegressionLocks(t *testing.T) {
	cases := []struct {
		pos, digit int
		short      bool // safe to run under -short
	}{
		{1000, 8, true},
		{10000, 7, true},
		{50000, 4, true},  // stdlib divider path (num is below the FFT gate post-truncation)
		{70000, 9, true},  // first to use the FFT divider (recip base case)
		{300000, 9, true}, // first to exercise recip's Newton recursion
		{100000, 4, false},
		{1000000, 5, false},
		{1000001, 1, false}, // famous millionth decimal digit
	}
	for _, c := range cases {
		t.Run(fmt.Sprintf("digit-%d", c.pos), func(t *testing.T) {
			if testing.Short() && !c.short {
				t.Skipf("skipping %d under -short", c.pos)
			}
			if got := extractDigit(c.pos); got != c.digit {
				t.Errorf("extractDigit(%d) = %d, want %d", c.pos, got, c.digit)
			}
		})
	}
}

// TestPow10 checks the 5^n·2^n decomposition against direct exponentiation,
// including a size where the 5-chain squarings cross the FFT threshold.
func TestPow10(t *testing.T) {
	for _, n := range []int{0, 1, 7, 32, 1000, 200000} {
		want := new(big.Int).Exp(big.NewInt(10), big.NewInt(int64(n)), nil)
		if pow10(n).Cmp(want) != 0 {
			t.Fatalf("pow10(%d) != 10^%d", n, n)
		}
	}
}

func TestEdges(t *testing.T) {
	if got := piDecimalString(0); got != "3" {
		t.Errorf("piDecimalString(0) = %q, want %q", got, "3")
	}
	if got := piDecimalString(1); got != "31" {
		t.Errorf("piDecimalString(1) = %q, want %q", got, "31")
	}
	if got := extractDigit(1); got != 3 {
		t.Errorf("extractDigit(1) = %d, want 3", got)
	}
	// digitPos below 1 is clamped to 1 by main; extractDigit(1) covers it.
}

// TestGuardBoundary asserts every digit up to the requested position matches
// the reference — the digits nearest the guard boundary are the first to rot
// if the precision/guard policy is cut too far.
func TestGuardBoundary(t *testing.T) {
	for _, d := range []int{100, 250, 500, 900, 1000} {
		s := piDecimalString(d)
		want := "3" + piRef[2:2+d]
		if s != want {
			t.Fatalf("guard boundary failure at d=%d", d)
		}
	}
}

// TestGuardStability is a reference-free guard against guardDigits being cut too
// small: ⌊π·10^d⌋ must be identical no matter how many guard digits we add. The
// digit at decimal place 4037 rots for guard < ~18, so this fails loudly if the
// guard regresses. (Position 4038 in the 1-based scheme = decimal place 4037.)
func TestGuardStability(t *testing.T) {
	for _, d := range []int{4038, 8000, 20000} {
		got := piFloorGuard(d, guardDigits, nil)
		ref := piFloorGuard(d, guardDigits+64, nil)
		if got.Cmp(ref) != 0 {
			t.Fatalf("guardDigits=%d insufficient at d=%d (differs from larger guard)", guardDigits, d)
		}
	}
}
