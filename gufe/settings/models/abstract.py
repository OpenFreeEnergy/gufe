from pydantic import BaseModel, Extra

class SettingsBaseModel(BaseModel):
    """
    Settings and modifications we want for all settings classes.
    """

    class Config:
        extra = Extra.forbid
        # Immutability in python is never strict.
        # If developers are determined/stupid
        # they can always modify a so-called "immutable" object.
        # https://pydantic-docs.helpmanual.io/usage/models/#faux-immutability
        allow_mutation = False
        arbitrary_types_allowed = True
        smart_union = True
