# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from .base import ExternalStorage

from ..errors import (
    MissingExternalResourceError, ChangedExternalResourceError
)

class S3Storage(ExternalStorage):
    """File storage backend for AWS S3.

    """
    def __init__(self, session: "boto3.Sessin", bucket: str, prefix: str):
        """

        """
        self.session = session

        self.resource = self.session.resource('s3')
        self.bucket = self.resource.Bucket(bucket)

        self.prefix = prefix

    def iter_contents(self, prefix=""):
        """Iterate over the labels in this storage.

        Parameters
        ----------
        prefix : str
            Only iterate over paths that start with the given prefix.

        Returns
        -------
        Iterator[str] :
            Contents of this storage, which may include items without
            metadata.
        """
        raise NotImplementedError()

    def _store_bytes(self, location, byte_data):
        """
        For implementers: This should be blocking, even if the storage
        backend allows asynchronous storage.
        """
        key = self.prefix + location

        self.bucket.put_object(Key=key, Body=byte_data)

    def _store_path(self, location, path):
        """
        For implementers: This should be blocking, even if the storage
        backend allows asynchronous storage.
        """
        """
        For implementers: This should be blocking, even if the storage
        backend allows asynchronous storage.
        """
        key = self.prefix + location

        with open(path, 'rb') as f:
            self.bucket.upload_fileobj(f, key)

    def _exists(self, location) -> bool:
        from botocore.exceptions import ClientError

        key = self.prefix + location

        # we do a metadata load as our existence check
        # appears to be most recommended approach
        try:
            self.bucket.Object(key).load()
            return True
        except ClientError:
            return False

    def _delete(self, location):
        key = self.prefix + location

        if self._exists(location):
            self.bucket.Object(key).delete()
        else:
            raise MissingExternalResourceError(
                f"Unable to delete '{str(key)}': Object does not exist"
            )

    def _get_filename(self, location):
        key = self.prefix + location

        object = self.bucket.Object(key)

        url = object.meta.client.generate_presigned_url(
                'get_object',
                ExpiresIn=0,
                Params={'Bucket': self.bucket.name, 'Key': object.key})

        # drop query params from url
        url = url.split('?')[0]

        return url

    def _load_stream(self, location):
        key = self.prefix + location

        try:
            return self.bucket.Object(key).get()['Body']
        except self.resource.meta.client.exceptions.NoSuchKey as e:
            raise MissingExternalResourceError(str(e))
