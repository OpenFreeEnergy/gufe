
class ProtocolUnitKey(str):
    def __repr__(self):
        return f"<ProtocolUnitKey('{str(self)}')>"


class ProtocolUnitMixin:

    @staticmethod
    def _keyencode_dependencies(inputs, targetcls):
        ninputs = dict()
        for key, value in inputs.items():
            if isinstance(value, dict):
                nvalue = dict()
                for k, v in value.items():
                    if isinstance(v, targetcls):
                        nvalue[k] = v.key
                    else:
                        nvalue[k] = v
            elif isinstance(value, list):
                nvalue = list()
                for i in value:
                    if isinstance(i, targetcls):
                        nvalue.append(i.key)
                    else:
                        nvalue.append(i)
            elif isinstance(value, targetcls):
                nvalue = value.key
            else:
                nvalue = value

            ninputs[key] = nvalue

        return ninputs
