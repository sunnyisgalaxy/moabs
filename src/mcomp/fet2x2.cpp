#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iterator>
#include <iostream>
#include <fstream>

#include "types.h"
#include "fisher_exact_test.h"

// Fisher 2x2 and 2x3 exact test command line utility
// Copyright (C) 2013  Christopher Chang  chrchang@alumni.caltech.edu

#define SMALLISH_EPSILON 0.00000000003
#define SMALL_EPSILON 0.0000000000001

// This helps us avoid premature floating point overflow.
#define EXACT_TEST_BIAS 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125

double fisher22(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t midp) {
  // Basic 2x2 Fisher exact test p-value calculation.
  double tprob = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
  double cur_prob = tprob;
  double cprob = 0;
  int32_t tie_ct = 1;
  uint32_t uii;
  double cur11;
  double cur12;
  double cur21;
  double cur22;
  double preaddp;
  // Ensure we are left of the distribution center, m11 <= m22, and m12 <= m21.
  if (m12 > m21) {
    uii = m12;
    m12 = m21;
    m21 = uii;
  }
  if (m11 > m22) {
    uii = m11;
    m11 = m22;
    m22 = uii;
  }
  if ((((uint64_t)m11) * m22) > (((uint64_t)m12) * m21)) {
    uii = m11;
    m11 = m12;
    m12 = uii;
    uii = m21;
    m21 = m22;
    m22 = uii;
  }
  cur11 = m11;
  cur12 = m12;
  cur21 = m21;
  cur22 = m22;
  while (cur12 > 0.5) {
    cur11 += 1;
    cur22 += 1;
    cur_prob *= (cur12 * cur21) / (cur11 * cur22);
    cur12 -= 1;
    cur21 -= 1;
    if (cur_prob < EXACT_TEST_BIAS) {
      if (cur_prob > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
        tie_ct++;
      }
      tprob += cur_prob;
      break;
    }
    cprob += cur_prob;
    if (cprob == INFINITY) {
      return 0;
    }
  }
  if ((cprob == 0) && (!midp)) {
    return 1;
  }
  if (cur12 > 0.5) {
    do {
      cur11 += 1;
      cur22 += 1;
      cur_prob *= (cur12 * cur21) / (cur11 * cur22);
      cur12 -= 1;
      cur21 -= 1;
      preaddp = tprob;
      tprob += cur_prob;
      if (tprob <= preaddp) {
	break;
      }
    } while (cur12 > 0.5);
  }
  if (m11) {
    cur11 = m11;
    cur12 = m12;
    cur21 = m21;
    cur22 = m22;
    cur_prob = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
    do {
      cur12 += 1;
      cur21 += 1;
      cur_prob *= (cur11 * cur22) / (cur12 * cur21);
      cur11 -= 1;
      cur22 -= 1;
      preaddp = tprob;
      tprob += cur_prob;
      if (tprob <= preaddp) {
        if (!midp) {
	  return preaddp / (cprob + preaddp);
        } else {
          return (preaddp - ((1 - SMALL_EPSILON) * EXACT_TEST_BIAS * 0.5) * tie_ct) / (cprob + preaddp);
        }
      }
    } while (cur11 > 0.5);
  }
  if (!midp) {
    return tprob / (cprob + tprob);
  } else {
    return (tprob - ((1 - SMALL_EPSILON) * EXACT_TEST_BIAS * 0.5) * tie_ct) / (cprob + tprob);
  }
}

double fisher22_1sided(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t m11_is_greater_alt, uint32_t midp) {
  double cur_prob = EXACT_TEST_BIAS;
  double left_prob = cur_prob;
  double right_prob = 0;
  uint32_t uii;
  double cur11;
  double cur12;
  double cur21;
  double cur22;
  double preaddp;
  // Ensure m11 <= m22 and m12 <= m21.
  if (m12 > m21) {
    uii = m12;
    m12 = m21;
    m21 = uii;
  }
  if (m11 > m22) {
    uii = m11;
    m11 = m22;
    m22 = uii;
  }
  // Flipping m11<->m12 and m21<->m22 also flips the direction of the
  // alternative hypothesis.  So we flip on m11-is-greater alternative
  // hypothesis here to allow the rest of the code to assume m11-is-less.
  if (m11_is_greater_alt) {
    uii = m11;
    m11 = m12;
    m12 = uii;
    uii = m21;
    m21 = m22;
    m22 = uii;
  }
  cur11 = m11;
  cur12 = m12;
  cur21 = m21;
  cur22 = m22;
  if ((((uint64_t)m11) * m22) >= (((uint64_t)m12) * m21)) {
    // starting right of (or at) center, p > 0.5
    // 1. left_prob = sum leftward to precision limit
    // 2. total_prob := left_prob
    // 3. total_prob += sum rightward to total_prob precision limit
    // return left_prob / total_prob
    while (cur11 > 0.5) {
      cur12 += 1;
      cur21 += 1;
      cur_prob *= (cur11 * cur22) / (cur12 * cur21);
      cur11 -= 1;
      cur22 -= 1;
      preaddp = left_prob;
      left_prob += cur_prob;
      if (left_prob <= preaddp) {
	break;
      }
      if (left_prob >= 1.0) {
	// Probability mass of our starting table was represented as 2^{-83},
	// so this would mean the left probability mass partial sum is greater
	// than 2^83 times that.  In which case the final p-value will
        // be indistinguishable from 1 at 53-bit precision if our input just
	// had 32-bit integers.  (Yes, the constant can be reduced.)
	return 1;
      }
    }
    cur11 = m11;
    cur12 = m12;
    cur21 = m21;
    cur22 = m22;
    cur_prob = EXACT_TEST_BIAS;
    right_prob = left_prob; // actually total_prob
    while (cur12 > 0.5) {
      cur11 += 1;
      cur22 += 1;
      cur_prob *= (cur12 * cur21) / (cur11 * cur22);
      cur12 -= 1;
      cur21 -= 1;
      preaddp = right_prob;
      right_prob += cur_prob;
      if (right_prob <= preaddp) {
	break;
      }
    }
    if (!midp) {
      return left_prob / right_prob;
    } else {
      return (left_prob - EXACT_TEST_BIAS * 0.5) / right_prob;
    }
  } else {
    // starting left of center, p could be small
    // 1. right_prob = sum rightward to precision limit
    // 2. left_prob = sum leftward to left_prob precision limit
    // return left_prob / (left_prob + right_prob)
    while (cur12 > 0.5) {
      cur11 += 1;
      cur22 += 1;
      cur_prob *= (cur12 * cur21) / (cur11 * cur22);
      cur12 -= 1;
      cur21 -= 1;
      preaddp = right_prob;
      right_prob += cur_prob;
      if (right_prob == INFINITY) {
	return 0;
      }
      if (right_prob <= preaddp) {
	break;
      }
    }
    cur11 = m11;
    cur12 = m12;
    cur21 = m21;
    cur22 = m22;
    cur_prob = EXACT_TEST_BIAS;
    while (cur11 > 0.5) {
      cur12 += 1;
      cur21 += 1;
      cur_prob *= (cur11 * cur22) / (cur12 * cur21);
      cur11 -= 1;
      cur22 -= 1;
      preaddp = left_prob;
      left_prob += cur_prob;
      if (left_prob <= preaddp) {
	break;
      }
    }
    if (!midp) {
      return left_prob / (left_prob + right_prob);
    } else {
      return (left_prob - EXACT_TEST_BIAS * 0.5) / (left_prob + right_prob);
    }
  }
}

void fisher22_precomp_thresh(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t* m11_minp, uint32_t* m11_maxp, uint32_t* tiep) {
  // Treating m11 as the only variable, this returns the minimum and (maximum -
  // 1) values of m11 which are less extreme than the observed result.  If the
  // observed result is maximally common, the return values will both be zero.
  // Also, if there is a second value of m11 which results in the exact same
  // p-value as the original, *tiep is set to that (otherwise it's just set to
  // the original m11).
  double cur_prob = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
  double cur11 = ((int32_t)m11);
  double cur12 = ((int32_t)m12);
  double cur21 = ((int32_t)m21);
  double cur22 = ((int32_t)m22);
  double ratio = (cur11 * cur22) / ((cur12 + 1) * (cur21 + 1));

  *tiep = m11;
  // Is m11 greater than the p-maximizing value?
  if (ratio > (1 + SMALL_EPSILON)) {
    *m11_maxp = m11;
    cur12 += 1;
    cur21 += 1;
    cur_prob *= ratio;
    do {
      cur11 -= 1;
      cur22 -= 1;
      m11--;
      cur12 += 1;
      cur21 += 1;
      cur_prob *= (cur11 * cur22) / (cur12 * cur21);
    } while (cur_prob > EXACT_TEST_BIAS);
    *m11_minp = m11;
    if (cur_prob > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
      *tiep = m11 - 1;
    }
  } else if (ratio > (1 - SMALL_EPSILON)) {
    *m11_minp = 0;
    *m11_maxp = 0;
    *tiep = m11 - 1;
  } else {
    // Is it less?
    cur11 += 1;
    cur22 += 1;
    ratio = (cur12 * cur21) / (cur11 * cur22);
    if (ratio > (1 + SMALL_EPSILON)) {
      cur_prob *= ratio;
      *m11_minp = ++m11;
      do {
	cur12 -= 1;
	cur21 -= 1;
	m11++;
	cur11 += 1;
	cur22 += 1;
	cur_prob *= (cur12 * cur21) / (cur11 * cur22);
      } while (cur_prob > EXACT_TEST_BIAS);
      *m11_maxp = m11;
      if (cur_prob > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
	*tiep = m11;
      }
    } else {
      *m11_minp = 0;
      *m11_maxp = 0;
      if (ratio > (1 - SMALL_EPSILON)) {
	*tiep = m11 + 1;
      }
    }
  }
}

int32_t fisher23_tailsum(double* base_probp, double* saved12p, double* saved13p, double* saved22p, double* saved23p, double *totalp, uint32_t* tie_ctp, uint32_t right_side) {
  double total = 0;
  double cur_prob = *base_probp;
  double tmp12 = *saved12p;
  double tmp13 = *saved13p;
  double tmp22 = *saved22p;
  double tmp23 = *saved23p;
  double tmps12;
  double tmps13;
  double tmps22;
  double tmps23;
  double prev_prob;
  // identify beginning of tail
  if (right_side) {
    if (cur_prob > EXACT_TEST_BIAS) {
      prev_prob = tmp13 * tmp22;
      while (prev_prob > 0.5) {
	tmp12 += 1;
	tmp23 += 1;
	cur_prob *= prev_prob / (tmp12 * tmp23);
	tmp13 -= 1;
	tmp22 -= 1;
	if (cur_prob <= EXACT_TEST_BIAS) {
	  break;
	}
	prev_prob = tmp13 * tmp22;
      }
      *base_probp = cur_prob;
      tmps12 = tmp12;
      tmps13 = tmp13;
      tmps22 = tmp22;
      tmps23 = tmp23;
    } else {
      tmps12 = tmp12;
      tmps13 = tmp13;
      tmps22 = tmp22;
      tmps23 = tmp23;
      while (1) {
	prev_prob = cur_prob;
	tmp13 += 1;
	tmp22 += 1;
	cur_prob *= (tmp12 * tmp23) / (tmp13 * tmp22);
	if (cur_prob < prev_prob) {
	  return 1;
	}
	tmp12 -= 1;
	tmp23 -= 1;
        if (cur_prob > (1 - 2 * SMALLISH_EPSILON) * EXACT_TEST_BIAS) {
	  // throw in extra (1 - SMALL_EPSILON) multiplier to prevent rounding
	  // errors from causing this to keep going when the left-side test
	  // stopped
	  if (cur_prob > (1 - SMALL_EPSILON) * EXACT_TEST_BIAS) {
	    break;
	  }
          *tie_ctp += 1;
	}
	total += cur_prob;
      }
      prev_prob = cur_prob;
      cur_prob = *base_probp;
      *base_probp = prev_prob;
    }
  } else {
    if (cur_prob > EXACT_TEST_BIAS) {
      prev_prob = tmp12 * tmp23;
      while (prev_prob > 0.5) {
	tmp13 += 1;
	tmp22 += 1;
	cur_prob *= prev_prob / (tmp13 * tmp22);
	tmp12 -= 1;
	tmp23 -= 1;
	if (cur_prob <= EXACT_TEST_BIAS) {
	  break;
	}
	prev_prob = tmp12 * tmp23;
      }
      *base_probp = cur_prob;
      tmps12 = tmp12;
      tmps13 = tmp13;
      tmps22 = tmp22;
      tmps23 = tmp23;
    } else {
      tmps12 = tmp12;
      tmps13 = tmp13;
      tmps22 = tmp22;
      tmps23 = tmp23;
      while (1) {
	prev_prob = cur_prob;
	tmp12 += 1;
	tmp23 += 1;
	cur_prob *= (tmp13 * tmp22) / (tmp12 * tmp23);
	if (cur_prob < prev_prob) {
	  return 1;
	}
	tmp13 -= 1;
	tmp22 -= 1;
        if (cur_prob > (1 - 2 * SMALLISH_EPSILON) * EXACT_TEST_BIAS) {
	  if (cur_prob > EXACT_TEST_BIAS) {
	    break;
	  }
	  *tie_ctp += 1;
	}
	total += cur_prob;
      }
      prev_prob = cur_prob;
      cur_prob = *base_probp;
      *base_probp = prev_prob;
    }
  }
  *saved12p = tmp12;
  *saved13p = tmp13;
  *saved22p = tmp22;
  *saved23p = tmp23;
  if (cur_prob > (1 - 2 * SMALLISH_EPSILON) * EXACT_TEST_BIAS) {
    if (cur_prob > EXACT_TEST_BIAS) {
      // even most extreme table on this side is too probable
      *totalp = 0;
      return 0;
    }
    *tie_ctp += 1;
  }
  // sum tail to floating point precision limit
  if (right_side) {
    prev_prob = total;
    total += cur_prob;
    while (total > prev_prob) {
      tmps12 += 1;
      tmps23 += 1;
      cur_prob *= (tmps13 * tmps22) / (tmps12 * tmps23);
      tmps13 -= 1;
      tmps22 -= 1;
      prev_prob = total;
      total += cur_prob;
    }
  } else {
    prev_prob = total;
    total += cur_prob;
    while (total > prev_prob) {
      tmps13 += 1;
      tmps22 += 1;
      cur_prob *= (tmps12 * tmps23) / (tmps13 * tmps22);
      tmps12 -= 1;
      tmps23 -= 1;
      prev_prob = total;
      total += cur_prob;
    }
  }
  *totalp = total;
  return 0;
}

double fisher23(uint32_t m11, uint32_t m12, uint32_t m13, uint32_t m21, uint32_t m22, uint32_t m23, uint32_t midp) {
  // 2x3 Fisher-Freeman-Halton exact test p-value calculation.
  // The number of tables involved here is still small enough that the network
  // algorithm (and the improved variants thereof that I've seen) are
  // suboptimal; a 2-dimensional version of the SNPHWE2 strategy has higher
  // performance.
  // 2x4, 2x5, and 3x3 should also be practical with this method, but beyond
  // that I doubt it's worth the trouble.
  // Complexity of approach is O(n^{df/2}), where n is number of observations.
  double cur_prob = (1 - SMALLISH_EPSILON) * EXACT_TEST_BIAS;
  double tprob = cur_prob;
  double cprob = 0;
  double dyy = 0;
  uint32_t tie_ct = 1;
  uint32_t dir = 0; // 0 = forwards, 1 = backwards
  double base_probl;
  double base_probr;
  double orig_base_probl;
  double orig_base_probr;
  double orig_row_prob;
  double row_prob;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double cur11;
  double cur21;
  double savedl12;
  double savedl13;
  double savedl22;
  double savedl23;
  double savedr12;
  double savedr13;
  double savedr22;
  double savedr23;
  double orig_savedl12;
  double orig_savedl13;
  double orig_savedl22;
  double orig_savedl23;
  double orig_savedr12;
  double orig_savedr13;
  double orig_savedr22;
  double orig_savedr23;
  double tmp12;
  double tmp13;
  double tmp22;
  double tmp23;
  double dxx;
  double preaddp;
  // Ensure m11 + m21 <= m12 + m22 <= m13 + m23.
  uii = m11 + m21;
  ujj = m12 + m22;
  if (uii > ujj) {
    ukk = m11;
    m11 = m12;
    m12 = ukk;
    ukk = m21;
    m21 = m22;
    m22 = ukk;
    ukk = uii;
    uii = ujj;
    ujj = ukk;
  }
  ukk = m13 + m23;
  if (ujj > ukk) {
    ujj = ukk;
    ukk = m12;
    m12 = m13;
    m13 = ukk;
    ukk = m22;
    m22 = m23;
    m23 = ukk;
  }
  if (uii > ujj) {
    ukk = m11;
    m11 = m12;
    m12 = ukk;
    ukk = m21;
    m21 = m22;
    m22 = ukk;
  }
  // Ensure majority of probability mass is in front of m11.
  if ((((uint64_t)m11) * (m22 + m23)) > (((uint64_t)m21) * (m12 + m13))) {
    ukk = m11;
    m11 = m21;
    m21 = ukk;
    ukk = m12;
    m12 = m22;
    m22 = ukk;
    ukk = m13;
    m13 = m23;
    m23 = ukk;
  }
  if ((((uint64_t)m12) * m23) > (((uint64_t)m13) * m22)) {
    base_probr = cur_prob;
    savedr12 = m12;
    savedr13 = m13;
    savedr22 = m22;
    savedr23 = m23;
    tmp12 = savedr12;
    tmp13 = savedr13;
    tmp22 = savedr22;
    tmp23 = savedr23;
    // m12 and m23 must be nonzero
    dxx = tmp12 * tmp23;
    do {
      tmp13 += 1;
      tmp22 += 1;
      cur_prob *= dxx / (tmp13 * tmp22);
      tmp12 -= 1;
      tmp23 -= 1;
      if (cur_prob <= EXACT_TEST_BIAS) {
        if (cur_prob > (1 - 2 * SMALLISH_EPSILON) * EXACT_TEST_BIAS) {
          tie_ct++;
	}
	tprob += cur_prob;
	break;
      }
      cprob += cur_prob;
      if (cprob == INFINITY) {
	return 0;
      }
      dxx = tmp12 * tmp23;
      // must enforce tmp12 >= 0 and tmp23 >= 0 since we're saving these
    } while (dxx > 0.5);
    savedl12 = tmp12;
    savedl13 = tmp13;
    savedl22 = tmp22;
    savedl23 = tmp23;
    base_probl = cur_prob;
    do {
      tmp13 += 1;
      tmp22 += 1;
      cur_prob *= (tmp12 * tmp23) / (tmp13 * tmp22);
      tmp12 -= 1;
      tmp23 -= 1;
      preaddp = tprob;
      tprob += cur_prob;
    } while (tprob > preaddp);
    tmp12 = savedr12;
    tmp13 = savedr13;
    tmp22 = savedr22;
    tmp23 = savedr23;
    cur_prob = base_probr;
    do {
      tmp12 += 1;
      tmp23 += 1;
      cur_prob *= (tmp13 * tmp22) / (tmp12 * tmp23);
      tmp13 -= 1;
      tmp22 -= 1;
      preaddp = tprob;
      tprob += cur_prob;
    } while (tprob > preaddp);
  } else {
    base_probl = cur_prob;
    savedl12 = m12;
    savedl13 = m13;
    savedl22 = m22;
    savedl23 = m23;
    if (!((((uint64_t)m12) * m23) + (((uint64_t)m13) * m22))) {
      base_probr = cur_prob;
      savedr12 = savedl12;
      savedr13 = savedl13;
      savedr22 = savedl22;
      savedr23 = savedl23;
    } else {
      tmp12 = savedl12;
      tmp13 = savedl13;
      tmp22 = savedl22;
      tmp23 = savedl23;
      dxx = tmp13 * tmp22;
      do {
	tmp12 += 1;
	tmp23 += 1;
	cur_prob *= dxx / (tmp12 * tmp23);
	tmp13 -= 1;
	tmp22 -= 1;
	if (cur_prob <= EXACT_TEST_BIAS) {
          if (cur_prob > (1 - 2 * SMALLISH_EPSILON) * EXACT_TEST_BIAS) {
	    tie_ct++;
	  }
	  tprob += cur_prob;
	  break;
	}
	cprob += cur_prob;
	if (cprob == INFINITY) {
	  return 0;
	}
	dxx = tmp13 * tmp22;
      } while (dxx > 0.5);
      savedr12 = tmp12;
      savedr13 = tmp13;
      savedr22 = tmp22;
      savedr23 = tmp23;
      base_probr = cur_prob;
      do {
	tmp12 += 1;
	tmp23 += 1;
	cur_prob *= (tmp13 * tmp22) / (tmp12 * tmp23);
	tmp13 -= 1;
	tmp22 -= 1;
	preaddp = tprob;
	tprob += cur_prob;
      } while (tprob > preaddp);
      tmp12 = savedl12;
      tmp13 = savedl13;
      tmp22 = savedl22;
      tmp23 = savedl23;
      cur_prob = base_probl;
      do {
	tmp13 += 1;
	tmp22 += 1;
	cur_prob *= (tmp12 * tmp23) / (tmp13 * tmp22);
	tmp12 -= 1;
	tmp23 -= 1;
	preaddp = tprob;
	tprob += cur_prob;
      } while (tprob > preaddp);
    }
  }
  row_prob = tprob + cprob;
  orig_base_probl = base_probl;
  orig_base_probr = base_probr;
  orig_row_prob = row_prob;
  orig_savedl12 = savedl12;
  orig_savedl13 = savedl13;
  orig_savedl22 = savedl22;
  orig_savedl23 = savedl23;
  orig_savedr12 = savedr12;
  orig_savedr13 = savedr13;
  orig_savedr22 = savedr22;
  orig_savedr23 = savedr23;
  for (; dir < 2; dir++) {
    cur11 = m11;
    cur21 = m21;
    if (dir) {
      base_probl = orig_base_probl;
      base_probr = orig_base_probr;
      row_prob = orig_row_prob;
      savedl12 = orig_savedl12;
      savedl13 = orig_savedl13;
      savedl22 = orig_savedl22;
      savedl23 = orig_savedl23;
      savedr12 = orig_savedr12;
      savedr13 = orig_savedr13;
      savedr22 = orig_savedr22;
      savedr23 = orig_savedr23;
      ukk = m11;
      if (ukk > m22 + m23) {
	ukk = m22 + m23;
      }
    } else {
      ukk = m21;
      if (ukk > m12 + m13) {
	ukk = m12 + m13;
      }
    }
    ukk++;
    while (--ukk) {
      if (dir) {
	cur21 += 1;
	if (savedl23) {
	  savedl13 += 1;
	  row_prob *= (cur11 * (savedl22 + savedl23)) / (cur21 * (savedl12 + savedl13));
	  base_probl *= (cur11 * savedl23) / (cur21 * savedl13);
	  savedl23 -= 1;
	} else {
	  savedl12 += 1;
	  row_prob *= (cur11 * (savedl22 + savedl23)) / (cur21 * (savedl12 + savedl13));
	  base_probl *= (cur11 * savedl22) / (cur21 * savedl12);
	  savedl22 -= 1;
	}
	cur11 -= 1;
      } else {
	cur11 += 1;
	if (savedl12) {
	  savedl22 += 1;
	  row_prob *= (cur21 * (savedl12 + savedl13)) / (cur11 * (savedl22 + savedl23));
	  base_probl *= (cur21 * savedl12) / (cur11 * savedl22);
	  savedl12 -= 1;
	} else {
	  savedl23 += 1;
	  row_prob *= (cur21 * (savedl12 + savedl13)) / (cur11 * (savedl22 + savedl23));
	  base_probl *= (cur21 * savedl13) / (cur11 * savedl23);
	  savedl13 -= 1;
	}
	cur21 -= 1;
      }
      if (fisher23_tailsum(&base_probl, &savedl12, &savedl13, &savedl22, &savedl23, &dxx, &tie_ct, 0)) {
	break;
      }
      tprob += dxx;
      if (dir) {
	if (savedr22) {
	  savedr12 += 1;
	  base_probr *= ((cur11 + 1) * savedr22) / (cur21 * savedr12);
	  savedr22 -= 1;
	} else {
	  savedr13 += 1;
	  base_probr *= ((cur11 + 1) * savedr23) / (cur21 * savedr13);
	  savedr23 -= 1;
	}
      } else {
	if (savedr13) {
	  savedr23 += 1;
	  base_probr *= ((cur21 + 1) * savedr13) / (cur11 * savedr23);
	  savedr13 -= 1;
	} else {
	  savedr22 += 1;
	  base_probr *= ((cur21 + 1) * savedr12) / (cur11 * savedr22);
	  savedr12 -= 1;
	}
      }
      fisher23_tailsum(&base_probr, &savedr12, &savedr13, &savedr22, &savedr23, &dyy, &tie_ct, 1);
      tprob += dyy;
      cprob += row_prob - dxx - dyy;
      if (cprob == INFINITY) {
	return 0;
      }
    }
    if (!ukk) {
      continue;
    }
    savedl12 += savedl13;
    savedl22 += savedl23;
    if (dir) {
      while (1) {
	preaddp = tprob;
	tprob += row_prob;
	if (tprob <= preaddp) {
	  break;
	}
	cur21 += 1;
	savedl12 += 1;
	row_prob *= (cur11 * savedl22) / (cur21 * savedl12);
	cur11 -= 1;
	savedl22 -= 1;
      }
    } else {
      while (1) {
	preaddp = tprob;
	tprob += row_prob;
	if (tprob <= preaddp) {
	  break;
	}
	cur11 += 1;
	savedl22 += 1;
	row_prob *= (cur21 * savedl12) / (cur11 * savedl22);
	cur21 -= 1;
	savedl12 -= 1;
      }
    }
  }
  if (!midp) {
    return tprob / (tprob + cprob);
  } else {
    return (tprob - ((1 - SMALLISH_EPSILON) * EXACT_TEST_BIAS * 0.5) * ((int32_t)tie_ct)) / (tprob + cprob);
  }
}

double fet2x2(int m11, int m12, int m21, int m22){
	uint32_t midp = 0; //no midp correction
	return fisher22(m11, m12, m21, m22, midp);
}


int fet_22( double & pvalue, int table[], int nrow, int ncol){
	if(nrow ==2 && ncol ==2){
		pvalue = fet2x2(table[0],table[1],table[2],table[3]);
	} else {
		fet_1k( pvalue, table, nrow, ncol );
	}
	return 0;
}
