
import logging
logger = logging.getLogger('iSDM.environment')
logger.setLevel(logging.DEBUG)

class EnvironmentalLayer(object):
    def __init__(self, **kwargs):
        # you want to be able to agregate at a different resolution
        # and back/forth, right?
        self.resolution = kwargs['resolution']

    def scale_resolution(degrees):
        pass

    def load_data(self, file_path=None):
        pass


class ClimateLayer(EnvironmentalLayer):
    pass

class LandCoverlayer(EnvironmentalLayer):
    pass

class LandUseLayer(EnvironmentalLayer):
    pass


