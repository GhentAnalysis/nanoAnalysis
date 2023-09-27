##########################
# Custom histogram class #
##########################
# Note:
#   This is just a wrapper around the builtin uproot class,
#   but with modifiable properties

class Histogram(object):

  def __init__(self, name=None, title=None, 
               edges=None, values=None, errors=None):
    self.name = name
    self.title = title
    self.edges = edges
    self.values = values
    self.errors = errors

  @classmethod
  def from_uproot(cls, hist):
    return cls(name=hist.name, title=hist.title,
               edges=hist.axis().edges(), values=hist.values(flow=True),
               errors=hist.errors(flow=True))

  def edges(self):
    return self.edges[:]

  def values(self):
    return self.values[:]

  def errors(self):
    return self.errors[:]
