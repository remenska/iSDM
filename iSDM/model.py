import os
from enum import Enum

import logging
logger = logging.getLogger('iSDM.model')
logger.setLevel(logging.DEBUG)


class Algorithm(Enum):
	GAM = 1
	GLM = 2
	MAXENT = 3


class Evaluation(Enum):
	ROC = 1
	KAPPA = 2
	TSS = 3


class Model(object):

	algorithm = None
	method = None
  species = None

	def __init__(self, algorithm = Algorithm.GAM):
		self.algorithm = algorithm

	def cross_validate(self, percentage, random_seed): 
		pass

	def evaluate_performance(self, method = Evaluation.ROC, **kwargs):
		self.method = method

	def fit(self):
		pass



