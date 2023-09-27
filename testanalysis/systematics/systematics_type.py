################################
# Get the type of a systematic #
################################

# Note:
#   The keys in this dict should be chosen correctly.
#   For weight systematics, the key must be the name of
#   one of the individual reweighters in the total reweighter being used
#   (see e.g. nanoAnalysis/reweighting/implementation/run2ulreweighter.py).
#   For selection systematics, the key must be in sync with systematics_tools.py.
#   For other systematics (currently not yet implemented),
#   the key must be recognized and correctly handled somewhere in eventloop.py

# Note:
#  They values in this dict should also be chosen correctly.
#  For weight systematics, use "weight", in order to indicate
#  that the variations and values should be calculated using the reweighter.
#  For selection systematics, use "[object]selection".
#  (For now, only the "selection" part is used and the "[object]" is ignored;
#   instead the key of the systematic is directly checked in eventloop.py
#   to decide which objects to recalculate; this might change in the future).
#  For other systematics: to be implemented.

systematics_type = ({
  #'muonreco': 'weight', # does not work yet
  'muonid': 'weight',
  'electronreco': 'weight',
  'electronid': 'weight',
  'pileup': 'weight',
  'prefire': 'weight',
  'btagging': 'weight',
  'njets': 'weight',
  'nbjets': 'weight',
  'jec': 'jetselection',
  'jer': 'jetselection',
  'uncl': 'metselection',
})
