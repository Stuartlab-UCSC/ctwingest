"""
Functions for interacting with the UCSC Cell Atlas Workbench API.
"""
import requests
from bs4 import BeautifulSoup
from ctwingest.io import read_json
from requests_toolbelt.multipart import encoder
import os
from requests.exceptions import ConnectionError


def read_credentials(creds_path):
    return read_json(creds_path)


def upload(ctw_path, credentials, url_base="http://localhost:5555/", group=None):
    assert credentials['email'] is not None and credentials["password"] is not None
    print("Attempting login")
    login_url = url_base + "user/sign-in"
    session = requests.Session()

    # Have to grab the token from the login request in order to successfully login.
    r = session.post(login_url, data=credentials)
    soup = BeautifulSoup(r.text, "html.parser")
    token = soup.find('input', {'name': 'csrf_token'})['value']
    credentials.update({"csrf_token": token})

    r = session.post(login_url, data=credentials)
    print("Successful login...")

    # The worksheet name is encoded by the filename of the .ctw.tgz file.
    worksheet_name = os.path.basename(ctw_path.split(".ctw.")[0])

    # Construct url.
    upload_url = url_base + "user/worksheet/" + worksheet_name
    if group is not None:
        upload_url += "?group=" + group

    # Tag the file onto the request.
    form = encoder.MultipartEncoder({
        "documents": (os.path.basename(ctw_path), open(ctw_path, "rb"), "application/octet-stream"),
        "composite": "NONE",
    })
    headers = {"Prefer": "respond-async", "Content-Type": form.content_type}

    # Attempt to upload file.
    try:
        r = session.post(upload_url, headers=headers, data=form)
    # TODO: the server isn't returning error codes with the requests library.
    # Catching this error is a temporary fix. Could be flask, could be machine based.
    except ConnectionError:
        # This error seems to occur in the likely case that
        # credentials aren't valid. TODO: evaluate whether the server
        # should be less cryptic, why isn't the 403 coming back?
        raise ValueError(
        """Please check your credentials file for the proper email and password.
        If you have supplied a group name, double check and make sure the name is spelled correctly."""
        )

    if r.status_code != 200:
        print("server failed with code %d" % r.status_code)
