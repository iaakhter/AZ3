/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Lev Nachmanson
*/
#pragma once
#include "util/lp/lp_settings.h"
#include "util/lp/lar_constraints.h"
namespace lean {
struct implied_bound {
    mpq m_bound;
    unsigned m_j; // the column for which the bound has been found
    bool m_is_low_bound;
    bool m_coeff_before_j_is_pos;
    unsigned m_row_or_term_index;
    bool m_strict;
    
    lconstraint_kind kind() const {
        lconstraint_kind k = m_is_low_bound? GE : LE;
        if (m_strict)
            k = static_cast<lconstraint_kind>(k / 2);
        return k;
    }
    bool operator==(const implied_bound & o) const {
        return m_j == o.m_j && m_is_low_bound == o.m_is_low_bound && m_bound == o.m_bound &&
            m_coeff_before_j_is_pos == o.m_coeff_before_j_is_pos &&
            m_row_or_term_index == o.m_row_or_term_index && m_strict == o.m_strict;
    }
    implied_bound(){}
    implied_bound(const mpq & a,
                  unsigned j,
                  bool low_bound,
                  bool coeff_before_j_is_pos,
                  unsigned row_or_term_index,
                  bool strict):
        m_bound(a),
        m_j(j),
        m_is_low_bound(low_bound),
        m_coeff_before_j_is_pos(coeff_before_j_is_pos),
        m_row_or_term_index(row_or_term_index),
        m_strict(strict) {}
};
}
