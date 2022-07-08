class ProtocolUnitToken(str):
    def __repr__(self):
        return f"<ProtocolUnitToken('{str(self)}')>"

class ProtocolUnitMixin:

    @staticmethod
    def _tokenize_dependencies(inputs, targetcls):
        ninputs = dict()
        for key, value in inputs.items():
            if isinstance(value, dict):
                nvalue = dict()
                for k, v in value.items():
                    if isinstance(v, targetcls):
                        nvalue[k] = v.token
                    else:
                        nvalue[k] = v
            elif isinstance(value, list):
                nvalue = list()
                for i in value:
                    if isinstance(v, targetcls):
                        nvalue.append(v.token)
                    else:
                        nvalue.append(v)
            elif isinstance(value, targetcls):
                nvalue = value.token
            ninputs[key] = nvalue

        return ninputs
