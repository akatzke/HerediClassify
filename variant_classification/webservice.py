#!/usr/bin/env python3

import pathlib
import datetime
import argparse

import uvicorn
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from classify import classify
from _version import __version__

from fastapi import Request, status
from fastapi.responses import JSONResponse
import traceback

import pybedtools

app = FastAPI()
@app.exception_handler(Exception)
async def my_exception_handler(request: Request, exc: Exception):
    return JSONResponse(
        status_code=500,
        content={"message": ''.join(traceback.format_exc()).replace('\n', '')}
    )

class Input(BaseModel):
    config_path: str
    variant_json: str


class Result(BaseModel):
    result: str
    config_file: str
    scheme_name: str
    scheme_version: str
    date: str
    tool_version: str


class Variant(BaseModel):
    variant: str


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
    variant_str = input.variant_json
    config_path = pathlib.Path(input.config_path)
    if not config_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"The config path {input.config_path} does not exist.",
        )
    final_config, classification_result = classify(config_path, variant_str)
    date = datetime.date.today().isoformat()
    pybedtools.cleanup()
    return Result(
        result=classification_result,
        config_file=config_path.name,
        scheme_name=final_config["name"],
        scheme_version=final_config["version"],
        date=date,
        tool_version=__version__
    )


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
