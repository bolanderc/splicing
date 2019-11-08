#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Splices together two near-field pressure signatures.

This routine splices together two near-field pressure signatures to aid in the
use of lower-fidelity tools for sections of a geometry that produce good
results. This can be used in conjunction with sBOOM and PyLdB to calculate the
loudness of the resulting sonic booms.

Routine Listings
-----------------

See Also
--------

Notes
-------

References
-----------

Example
--------

"""

import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
from rapidboom.sboomwrapper import SboomWrapper
import pyldb.core as pyldb


class Splicing:
    def __init__(self, front_sig, rear_sig, x_c, blending, l_ref):
        self.x_c = x_c
        self.x_front = front_sig[:, 0]
        self.dp_p_front = front_sig[:, 1]
        self.n_front = len(self.dp_p_front)
        self.x0 = self.x_front[0]
        self.x_rear = rear_sig[:, 0]
        self.dp_p_rear = rear_sig[:, 1]
        self.n_rear = len(self.dp_p_rear)
        self.n_blend = blending[0]
        self.blend_type = blending[1]
        self.l_ref = l_ref
        self.dp_p_front[0] = 0.
        self.dp_p_rear[0] = 0.
        if self.blend_type is 'linear':
            self.blend_param = np.linspace(0., 1., num=self.n_blend)
        elif self.blend_type is 'quad':
            self.blend_param = np.linspace(0., 1., num=self.n_blend)
        elif self.blend_type is 'cosine':
            self.blend_param = np.linspace(0., np.pi/2., num=self.n_blend)
        else:
            print('Blending type not found')
            return

    def propagate_sig(self, case, nearfield_sig, R_over_L=3):
        REF_LENGTH = self.l_ref
        sboom = SboomWrapper('./', exe='./sboom_linux')
        sboom.set(signature=nearfield_sig,
                  mach_number=1.6,
                  altitude=51706.037,
                  propagation_start=R_over_L*REF_LENGTH*3.28084,
                  altitude_stop=0.,
                  output_format=0,
                  input_xdim=2,
                  propagation_points=40000,
                  padding_points=8000)

        results = sboom.run()
        g_sig = results["signal_0"]["ground_sig"]
        np.savetxt("ground_"+case, g_sig)

        self.noise_level = pyldb.perceivedloudness(g_sig[:, 0], g_sig[:, 1],
                                                   pad_rear=4)
        print(self.noise_level)
        return self.noise_level

    def splice_sigs(self, n_window=50, save_sig=[True, 'nearfield'],
                    plot=False, propagate=[False, 'nearfield', 3]):
#        self._crop_sigs(n_window)
        print('w')
        self._nondimensionalize_sigs()
        self._cut_and_blend()
        spliced_nf_xnd = np.concatenate((self.x_front_cutnd,
                                         self.x_rear_cutnd))
        spliced_nf_x = spliced_nf_xnd*self.l_ref + self.x0
        spliced_nf_dp = np.concatenate((self.dp_front_blend,
                                        self.dp_rear_cut))
        spliced_nfnd = np.array([spliced_nf_xnd, spliced_nf_dp]).T
        spliced_nf = np.array([spliced_nf_x, spliced_nf_dp]).T
        if plot:
            plt.plot(spliced_nfnd[:, 0], spliced_nfnd[:, 1])
            plt.title("Spliced Near-field Signature (Nondim)")
            plt.tight_layout()
            plt.show()
            plt.plot(spliced_nf[:, 0], spliced_nf[:, 1])
            plt.title("Spliced Near-field Signature")
            plt.tight_layout()
            plt.show()
        if save_sig[0]:
            np.savetxt(save_sig[1], spliced_nf)
            np.savetxt(save_sig[1] + '_nondim', spliced_nfnd)
        if propagate[0]:
            self.propagate_sig(propagate[1], spliced_nf, propagate[2])
        return spliced_nf

    def _nondimensionalize_sigs(self):
        self.x_front_nd = np.array((self.x_front - self.x0)/self.l_ref)
        self.x_rear_nd = np.array((self.x_rear - self.x_rear[0])/self.l_ref)
        self.dp_p_front_interp = interpolate.interp1d(self.x_front_nd,
                                                      self.dp_p_front,
                                                      kind='cubic')
        self.dp_p_rear_interp = interpolate.interp1d(self.x_rear_nd,
                                                     self.dp_p_rear,
                                                     kind='cubic')

    def _crop_sigs(self, n_window):
        p_start_f = self.x_front[self.dp_p_front > 1000.*self.dp_p_front[0]][0]
        p_start_r = self.x_rear[self.dp_p_rear > 1000.*self.dp_p_rear[0]][0]
        n_pad_f = len(self.x_front[self.x_front < p_start_f])
        n_pad_r = len(self.x_rear[self.x_rear < p_start_r])
        num_x_del = np.abs(n_pad_f - n_pad_r)
        i_range_del = np.arange(0, num_x_del + 1)
        if p_start_f < p_start_r:
            self.x_rear = np.delete(self.x_rear, i_range_del)
            self.dp_p_rear = np.delete(self.dp_p_rear, i_range_del)
            dx = self.x_rear[0] - self.x_front[0]
            self.x_rear -= dx
        else:
            self.x_front = np.delete(self.x_front, i_range_del)
            self.dp_p_front = np.delete(self.dp_p_front, i_range_del)
            dx = self.x_front[0] - self.x_rear[0]
            self.x_front -= dx
        window = np.cos(np.linspace(0., np.pi/2., n_window))
        self.dp_p_front[-n_window:] *= window
        if self.x_front[-1] < self.x_rear[-1]:
            dx = self.x_front[1] - self.x_front[0]
            x_pad = np.arange(self.x_front[-1] + dx, self.x_rear[-1] + dx, dx)
            dp_pad = np.zeros_like(x_pad)
            self.x_front = np.concatenate((self.x_front, x_pad))
            self.dp_p_front = np.concatenate((self.dp_p_front, dp_pad))
        self.x_front[-1] = self.x_rear[-1]
        self.x_front[0] = self.x_rear[0]

    def _cut_and_blend(self):
        front_cut_ind = np.where((self.x_rear_nd <= self.x_c))
        self.x_front_cutnd = self.x_rear_nd[front_cut_ind]
        self.dp_front_cut = self.dp_p_front_interp(self.x_front_cutnd)
        rear_cut_ind = np.where((self.x_rear_nd >= self.x_c))
        self.x_rear_cutnd = self.x_rear_nd[rear_cut_ind]
        self.dp_rear_cut = self.dp_p_rear_interp(self.x_rear_cutnd)
        self.dp_front_blend = np.copy(self.dp_front_cut)
        if self.blend_type is 'linear':
            blend_const = [(1. - self.blend_param), self.blend_param]
        elif self.blend_type is 'quad':
            blend_const = [(1. - self.blend_param**2), self.blend_param**2]
        elif self.blend_type is 'cosine':
            blend_const = [np.cos(self.blend_param), np.sin(self.blend_param)]
        self.dp_front_blend[-self.n_blend:] *= blend_const[0]
        self.dp_front_blend[-self.n_blend:] += self.dp_rear_cut[0]*blend_const[1]
