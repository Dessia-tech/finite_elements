import pydocstyle
import os
from glob import glob
import random
from datetime import date

file_list = filter(lambda z: not z.endswith("__init__.py"),
                   [y for x in os.walk('./aa_toulouse')
                    for y in glob(os.path.join(x[0], '*.py'))])

UNWATCHED_ERRORS = [
    # Do not watch these errors
    'D100', 'D104', 'D105', 'D107',
    'D200', 'D202', 'D203', 'D204', 'D205', 'D206', 'D210', 'D212',
    'D301', 'D302',
    'D400', 'D401', 'D402', 'D407', 'D408', 'D409',
    'D412', 'D415', 'D418'
]

MAX_ERROR_BY_TYPE = {
    # http://www.pydocstyle.org/en/stable/error_codes.html
    'D100': 0,
    'D101': 0,
    'D102': 0,
    'D103': 0,
    'D104': 0,
    'D105': 0,
    'D106': 0,
    'D107': 0,

    'D200': 0,
    'D201': 0,
    'D202': 0,
    'D203': 0,
    'D204': 0,
    'D206': 0,
    'D207': 0,
    'D208': 0,
    'D209': 0,
    'D210': 0,
    'D211': 0,
    'D212': 0,
    'D213': 0,
    'D214': 0,
    'D215': 0,

    'D300': 0,
    'D301': 0,
    'D302': 0,

    'D401': 0,
    'D402': 0,
    'D403': 0,
    'D404': 0,
    'D405': 0,
    'D406': 0,
    'D407': 0,
    'D408': 0,
    'D409': 0,
    'D410': 0,
    'D411': 0,
    'D412': 0,
    'D413': 0,
    'D414': 0,
    'D415': 0,
    'D416': 0,
    'D417': 0,
    'D418': 0,
}

error_detected = False
error_over_ratchet_limit = False
ratchet_limit = 9
effective_date = date(2022, 12, 20)
today = date.today()
weekly_decrease = 5
time_decrease = int((today - effective_date).days / 7. * weekly_decrease)


code_to_errors = {}
for error in pydocstyle.check(file_list, ignore=UNWATCHED_ERRORS):
    code_to_errors.setdefault(error.code, [])
    code_to_errors[error.code].append(error)

code_to_number = {code: len(errors) for code, errors in code_to_errors.items()}

for error_code, number_errors in code_to_number.items():
    if error_code not in UNWATCHED_ERRORS:
        max_errors = max(MAX_ERROR_BY_TYPE.get(error_code, 0) - time_decrease, 0)

        if number_errors > max_errors:
            error_detected = True
            print(f'\nFix some {error_code} errors: {number_errors}/{max_errors}')

            errors = code_to_errors[error_code]
            errors_to_show = sorted(random.sample(errors, min(30, len(errors))),
                                    key=lambda m: (m.filename, m.line))
            for error in errors_to_show:
                print(f'{error.filename} line {error.line}: {error.message}')
        elif max_errors - ratchet_limit <= number_errors < max_errors:
            print(
                f'\nYou can lower number of {error_code} to {number_errors + time_decrease} (actual {max_errors + time_decrease})')
        elif number_errors < max_errors - ratchet_limit:
            error_over_ratchet_limit = True
            print(
                f'\nYou MUST lower number of {error_code} to {number_errors + time_decrease} (actual {max_errors + time_decrease})')

if error_detected:
    raise RuntimeError('Too many errors\nRun pydocstyle aa_toulouse to get the errors')

if error_over_ratchet_limit:
    raise RuntimeError(
        'Please lower the error limits in code_pydocstyle.py MAX_ERROR_BY_TYPE according to warnings above')
