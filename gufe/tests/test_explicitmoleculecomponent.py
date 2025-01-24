import pickle

class ExplicitMoleculeComponentMixin:

    def test_pickle(self, instance):

        pickled = pickle.dumps(instance)
        unpickled = pickle.loads(pickled)

        assert unpickled == instance

        # it's currently the case that the flyweight pattern isn't respected by
        # pickling/unpickling
        assert not unpickled is instance
