import urllib.parse
import flask
import urllib
import os
import HorusAPI

load_page = HorusAPI.PluginPage(
    id="load_tables",
    name="Load tables",
    description="Load tables",
    html="index.html",
    hidden=True,
)


def load_html():

    path = flask.request.args.get("path")

    # The path was safely encoded in the uri
    path = urllib.parse.unquote(path)

    if path is None or not os.path.exists(path):
        return flask.jsonify({"error": "Path does not exist"}, 404)

    return flask.send_file(path)


load_page.addEndpoint(
    HorusAPI.PluginEndpoint(url="/load_table", methods=["GET"], function=load_html)
)
