from gufe import validation


def test_validation(alchem_network_benzene_variants):
    validation.chain_validate(alchem_network_benzene_variants)
