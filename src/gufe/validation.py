from gufe.tokenization import GufeTokenizable, KeyedChain

# def chain_validate(gt: GufeTokenizable) -> tuple[set[GufeTokenizable], set[GufeTokenizable, Exception], set[GufeTokenizable]]:
#     valid: set[GufeTokenizable] = set()
#     invalid: set[tuple[GufeTokenizable, Exception | None]] = set()
#     not_eligible: set[GufeTokenizable]= set()

#     for key, token in gt.to_keyed_chain():
#         if hasattr(token, "validate"):
#             is_valid = token.validate()

#     return (valid, invalid, not_eligible)


def chain_validate(gt: GufeTokenizable):
    for token in KeyedChain.from_gufe(gt).decode_subchains(func=lambda x: True):
        if hasattr(token, "is_valid"):
            token.is_valid()
            continue
        print(token)
    return None
