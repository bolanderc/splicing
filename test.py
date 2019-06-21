#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pyldb unit testing
"""

import pytest
import splicing
import numpy as np
import os.path


def test_splicing():
    front_test_sig_fname = os.path.join("misc", "25DNosePAN")
    rear_test_sig_fname = os.path.join("misc", "axie_panair_nf.txt")
    front_sig = np.genfromtxt(front_test_sig_fname)
    rear_sig = np.genfromtxt(rear_test_sig_fname)
    splice = splicing.Splicing(front_sig,
                               rear_sig,
                               0.4377, [10, 'quad'],
                               32.92)
    nf = splice.splice_sigs(save_sig=[False, '25D_splice'])
    sig_int = np.trapz(nf[:, 1], x=nf[:, 0])
    sig_int_test = 0.022264520953522857

    assert np.allclose(sig_int, sig_int_test, rtol=0.0, atol=10e-12)
