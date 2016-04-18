
class EnvironmentalLayer(object):
    def __init__(self, **kwargs):
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


