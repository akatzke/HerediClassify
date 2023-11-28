#!/usr/bin/env python3

import pathlib
import datetime
import argparse

import uvicorn
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from classify import classify
from _version import __version_

app = FastAPI()


class Input(BaseModel):
    config_path: str
    variant_json: str


class Result(BaseModel):
    result: str
    config: str
    date: str
    version: str == __version__


@app.get("_ping")
async def _ping():
    pass


summary = "Classify variant"


@app.post(
    "/classify_variant",
    summary=summary,
    description=f"""
          {summary}

          Parameters
          ---------
          variant: str
              Json string containing the variant information
          config_path: str
              Path to classification config


          Returns
          ---------
          result: str
              Json string containing the variant classification result
          config: str
              Name of classification configuration file
          date: str
              Date of automated classification
          version: str
              The version of the classification tool

          """,
)
async def classify_variant(input: Input) -> Result:
    """
    Execute classification of variant
    """
    variant_json = input.variant_json
    config_path = pathlib.Path(input.config_path)
    if not config_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"The config path {input.config_path} does not exist.",
        )
    classification_result = classify(config_path, variant_json)
    date = datetime.date.today().isoformat()
    return Result(classification_result, config_path.name, date)


def main():
    # define CLI arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--port", action="store", default=8080, help="Port to listen on", type=int
    )
    parser.add_argument(
        "--host", action="store", default="0.0.0.0", help="Hosts to listen on", type=str
    )

    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )

    # read passed CLI arguments
    args = parser.parse_args()

    # create and run the web service
    uvicorn.run(app, host=args.host, port=args.port, reload=False)


if __name__ == "__main__":
    main()
