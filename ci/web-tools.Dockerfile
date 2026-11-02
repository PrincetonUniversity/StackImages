# Tools for the web page: lint, type checks, unit + Chromium tests, Pyodide vendoring,
# site and snapshot builds. Python 3.13 matches Pyodide 0.29.3.
FROM python:3.13.15-slim-bookworm
ENV PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLCONFIGDIR=/tmp/matplotlib \
    PLAYWRIGHT_BROWSERS_PATH=/ms-playwright
COPY requirements-web.lock /tmp/requirements-web.lock
RUN pip install -r /tmp/requirements-web.lock \
 && python -m playwright install --with-deps chromium
WORKDIR /work
CMD ["python", "-m", "pytest", "-v", "-p", "no:cacheprovider", "tests/web"]
