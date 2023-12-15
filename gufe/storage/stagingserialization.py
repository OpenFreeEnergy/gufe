from gufe.tokenization import JSON_HANDLER
from gufe.custom_json import JSONCodec, JSONSerializerDeserializer
from .stagingregistry import StagingPath


class StagingPathSerialization:
    """Class for managing serialization of a :class:`.StagingPath`.

    Serialization of a :class:`.StagingPath` needs to strip the specific
    storage context (path to external files storage) because we should able
    to change that out-of-process (e.g., move the directory containing
    results) and still be able to deserialize correctly. This class is
    responsible for abstracting/injecting the storage context for a
    :class:`.StagingPath` is serialized/deserialized.
    """
    # TODO: this long comment should probably go somewhere where it will
    # show up in docs as well? Maybe just bump it into the class docstring?
    #
    # Serializing staging paths
    # -------------------------
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
    # 5. I have generated data in one backend, and I tranferred it to
    #    another backend. It needs to be readable from the other backend.
    #    (Use case: data is in long-term cloud storage that requires
    #    credentials, but I want to share some part of that data with
    #    someone else by transferring it to a disk.)
    # 6. I am interfacing with a package that adds serialization types to
    #    the gufe JSON_HANDLER via an external JSONCodec. Maybe, in the
    #    worst case, the external codec gets added *after* I've created my
    #    serialization object. I need to be able to serialize those custom
    #    types.
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
    # User stories 3--6 are handled by this
    # :class:`.StagingPathSerialization` class. Story 3 is handled by
    # allowing the appropriate context (in the form of a
    # :class:`.StorageManager`) to be injected into the deserialization
    # process. Story 4 can be handled by using more than one
    # :class:`.StagingPathSerialization` context (associated with different
    # :class:`.StorageManager` objects. Story 5 is handled by injecting
    # the appropriate context (and, in principle, is a variant of story 3.)
    # Story 6 is handled by doing a just-in-time generation of the
    # JSONSerializerDeserializer that we use for this class.
    def __init__(self, manager):
        self.manager = manager
        self.codec = JSONCodec(
            cls=StagingPath,
            to_dict=self.to_dict,
            from_dict=self.from_dict,
        )
        self.refresh_handler()

    def refresh_handler(self):
        """Ensure that the current handler includes all registered codecs"""
        codecs = JSON_HANDLER.codecs + [self.codec]
        self.json_handler = JSONSerializerDeserializer(codecs)

    @property
    def encoder(self):
        """
        JSONEncoder class to use when serializing a :class:`.StagingPath`
        """
        self.refresh_handler()
        return self.json_handler.encoder

    @property
    def decoder(self):
        """
        JSONdecoder class to use when deserializing a :class:`.StagingPath`
        """
        self.refresh_handler()
        return self.json_handler.decoder

    def to_dict(self, path: StagingPath):
        """
        Dict representation of a StagingPath, abstracting specific context.

        This provides a JSON-serializable representation of a StagingPath
        where the specific context of the StagingPath (the specific storage
        backend where it is located) is replaced by a generic representation
        of 'scratch', 'shared', or 'permanent', allowing a new specific
        context to be injected on deserialization.
        """
        # scratch, shared, permanent may form nested with progressively
        # smaller contexts, so the last of those it is in is where it should
        # be labelled. TODO: opportunity for performance improvement if
        # needed
        loc = None
        if path.label in self.manager.scratch_root.iterdir():
            # TODO: does this happen? we should only trigger this function
            # on a StagingPath, and anything in scratch will only be
            # pathlib.Path, right?
            loc = "scratch"
        if path.label in self.manager.shared_root.iter_contents():
            loc = "shared"
        if path.label in self.manager.permanent_root.iter_contents():
            loc = "permanent"

        if loc is None:
            raise RuntimeError(
                f"Unable to serialize {path}: it does not appear to be "
                "associated with storage managed by the context manager "
                f"{self.manager}."
            )

        return {
            ':container:': loc,
            ':label:': path.label,
        }

    def from_dict(self, dct: dict) -> StagingPath:
        """Recreate a StagingPath from its dict represnetation.

        This undoes the process from :method:`.to_dict`. It injects the
        storage context in ``self.storage_manager`` into the deserialized
        :class:`.StagingPath` instance.
        """
        staging = getattr(self.manager, f"{dct[':container:']}_staging")
        return staging / dct[':label:']
