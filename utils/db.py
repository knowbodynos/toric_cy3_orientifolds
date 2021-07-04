from pymongo import MongoClient
from pymongo.database import Database


MONGO_URI = "mongodb://frontend:password@129.10.135.232:27017/ToricCY"


class ToricCY(Database):
    def __init__(self, uri: str = MONGO_URI):
        client = MongoClient(uri)
        self.__dict__ = client.ToricCY.__dict__

    def find_invols(self, *args, **kwargs):
        invols = list(self.INVOL.find(*args, **kwargs))
        for x in invols:
            x.update(self.TRIANG.find_one({k: x[k] for k in ["POLYID", "GEOMN", "TRIANGN"]}))
            x.update(self.GEOM.find_one({k: x[k] for k in ["POLYID", "GEOMN"]}))
            x.update(self.POLY.find_one({"POLYID": x["POLYID"]}))
        return invols
