from gufe.tokenization import JSON_HANDLER
from gufe.custom_json import JSONCodec, JSONSerializerDeserializer
from .stagingdirectory import StagingPath


class StagingPathSerialization:
    # TODO: where should this go? I think maybe on the storage manager

    # Serializing staging paths
    # -------------------------
    #
    # Some important user stories to consider:
    #
    # 1. I am loading my results object, and I will want to use the
    #    associated files. This should be transparent, regardless of where
    #    the permanent storage is located.
    # 2. I am loading my results object, but I do not need the large stored
    #    files. I do not want to download them when they aren't needed.
    # 3. My permanent storage was a directory on my file system, but I have
    #    moved that directory (with use cases of (a) I moved the absolute
    #    path; (b) it is at a different relative path with respect to my
    #    pwd.
    # 4. I'm working with files from two different permanent storages. I
    #    need to be able to load from both in the same Python process.
    #
    # Outputs from a protocol may contain a :class:`.StagingPath`. Note that
    # a :class:`.StagingPath` is inherently not a complete description of
    # how to obtain the associated data: in almost every situation, there is
    # some additional context required. This can include the credentials to
    # an external server, or a base path that the file can be found relative
    # to (which may have changed if the user moves it.) Because of this, we
    # need to inject that context in to the deserialization.
    #
    # This object injects the relevant context, provided by the
    # :class:`.StorageManager`. It creates a JSONSerializerDeserializer
    # based on the one being used by gufe in this process, using all the
    # installed codecs plus an additional codec specific to this context.
    #
    # User stories 1 and 2 are handled by the nature of the
    # :class:`.StagingPath` object. The external file downloads as part of
    # the ``__fspath__`` method. This means that when using the ``open``
    # builtin, you will automatically download the file to a local staging
    # directory. However, the reference to the file can exist in the results
    # object without downloading the file.
    #
    # User stories 3 and 4 are handled by this
    # :class:`.StagingPathSerialization` class. Story 3 is handled by
    # allowing the appropriate context (in the form of a
    # :class:`.StorageManager`) to be injected into the deserialization
    # process. Story 4 can be handled by using more than one
    # :class:`.StagingPathSerialization` context (associated with different
    # :class:`.StorageManager` objects.

    def __init__(self, manager):
        self.manager = manager
        self.codec = JSONCodec(
            cls=StagingPath,
            to_dict=self.to_dict,
            from_dict=self.from_dict,
        )
        self.refresh_handler()

    def refresh_handler(self):
        codecs = JSON_HANDLER.codecs + [self.codec]
        self.json_handler = JSONSerializerDeserializer(codecs)

    @property
    def encoder(self):
        return self.json_handler.encoder

    @property
    def decoder(self):
        return self.json_handler.decoder

    def to_dict(self, path):
        # scratch, shared, permanent may form nested with progressively
        # smaller contexts, so the last of those it is in is where it should
        # be labelled. TODO: opportunity for performance improvement if
        # needed
        loc = None
        if path.label in self.manager.scratch_root.iterdir():
            loc = "scratch"
        if path.label in self.manager.shared_root.iter_contents():
            loc = "shared"
        if path.label in self.manager.permanent_root.iter_contents():
            loc = "permanent"

        return {
            ':container:': loc,
            ':label:': path.label,
        }

    def from_dict(self, dct):
        staging = getattr(self.manager, f"{dct[':container:']}_staging")
        return staging / dct[':label:']
