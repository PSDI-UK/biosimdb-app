#!/usr/bin/env python
"""
gunicorn WSGI server configuration.
"""

secure_scheme_headers = {
    'X-FORWARDED-PROTOCOL': 'http',
    'X-FORWARDED-PROTO': 'http',
    'X-FORWARDED-SSL': 'off'
}