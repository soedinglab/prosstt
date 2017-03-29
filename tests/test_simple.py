import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from prosstt import simulation as sim

import unittest
import numpy as np

class TestSimulationMethods(unittest.TestCase):

	def test_are_lengths_ok(self):
		# when Ms is empty it should return false
		Ms = None
		self.assertFalse(sim.are_lengths_ok(Ms))
		# max_crit is false
		Ms = [np.ones((50, 50)) for i in range(3)]
		Ms = Ms*500
		Ms[0][0][0] = 900
		self.assertFalse(sim.are_lengths_ok(Ms))
		# rel_crit is false
		Ms = [np.ones((50, 50)) for i in range(3)]
		Ms[0][0][0] = 700
		self.assertFalse(sim.are_lengths_ok(Ms))
