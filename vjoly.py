import bz2
import collections
import csv
import datetime
import gzip
import os
import re
import subprocess
import sys

import Bio.Alphabet
import Bio.Alphabet.IUPAC
import Bio.SeqIO
import Bio.SeqRecord
import tabulate


# CONSTANTS
STDIO_PATH = '-'

BAM_EXTS = ('.bam', )
FASTA_EXTS = ('.fa', '.fas', '.fasta', '.faa', '.fna')
FASTQ_EXTS = ('.fq', '.fastq')
GENBANK_EXTS = ('.gb', '.gbff', '.seq')
SEQFILE_EXTS = FASTA_EXTS + GENBANK_EXTS

GZIP_EXT = '.gz'
BGZIP_EXT = '.bgz'
BZIP2_EXT = '.bz2'
ZIP_EXTS = (GZIP_EXT, BGZIP_EXT, BZIP2_EXT)

FASTA_FORMAT = 'fasta'
GENBANK_FORMAT = 'genbank'
SEQFILE_FORMAT = FASTA_FORMAT
FILE_INDEXED = False

NT_ALPHA = Bio.Alphabet.NucleotideAlphabet
AA_ALPHA = Bio.Alphabet.ProteinAlphabet
DNA_ALPHA = Bio.Alphabet.IUPAC.IUPACUnambiguousDNA
DNA_ALPHA_EXT = Bio.Alphabet.IUPAC.IUPACAmbiguousDNA
PROT_ALPHA = Bio.Alphabet.IUPAC.IUPACProtein
PROT_ALPHA_EXT = Bio.Alphabet.IUPAC.ExtendedIUPACProtein
RNA_ALPHA = Bio.Alphabet.IUPAC.IUPACUnambiguousRNA
RNA_ALPHA_EXT = Bio.Alphabet.IUPAC.IUPACAmbiguousRNA


# MESSAGE/ERROR DISPLAY
def print_stderr(msg, date=False, time=False, program=None,
                 prog_width=20, indent=0, indent_width=4,
                 wrap=80, min_width=25, prefix=None, end='\n'):
    left = []
    if date or time:
        now = datetime.datetime.now()
        if date and time:
            left.append('{:%Y-%m-%d %H:%M:%S}'.format(now))
        elif date:
            left.append('{:%Y-%m-%d}'.format(now))
        else:
            left.append('{:%H:%M:%S}'.format(now))
    if program is not None:
        program = '{:<{w}.{w}}'.format(program, w=prog_width)
        left.append('[{}]'.format(program))

    left = '  '.join(left)
    w = indent * indent_width
    if left:
        w += len(left) + 2
    if wrap - w < min_width:
        left += '\n'
        indent += 1
        w = indent * indent_width
    elif left:
        left += '  '
    wrap -= w
    sep = '\n' + ' ' * w

    wrap_regex = '.{{1,{:d}}}(?:\\s+|$)'.format(wrap)
    if prefix:
        msg = '{}: {}'.format(prefix, msg)
    right = sep.join([line.strip() for line in re.findall(wrap_regex, msg)])
    right = ' ' * indent * indent_width + right
    sys.stderr.write(left + right + end)


def print_table(table, header=(), tablefmt='simple'):
    print('')
    print(tabulate.tabulate(table, headers=header, tablefmt=tablefmt))
    print('')


def multiple_choice(question, options, default=None, indent=0):
    options = [i.lower() for i in options]
    default = default.lower()
    options_str = []
    for option in options:
        if default is not None and option == default:
            options_str.append(option.upper())
        else:
            options_str.append(option.lower())
    options_str = '[{}]'.format('/'.join(options_str))
    print_stderr(question + ' ' + options_str, end=' ', indent=indent)
    while True:
        choice = input().lower()
        if not choice and default is not None:
            choice = default
        if choice in options:
            break
        valid_options = '[{}] and [{}]'.format(options[-2], options[-1])
        if len(options) > 2:
            first_options = ', '.join(['[' + i + ']' for i in options[:-2]])
            valid_options = first_options + ', ' + valid_options
        msg = 'Please choose between {}:'.format(valid_options)
        print_stderr(msg, end=' ', indent=indent)
    return choice


# FILE UTILS
def open_file(path, mode='r'):
    if mode not in ('r', 'w'):
        raise ValueError("Argument 'mode' must be 'r' or 'w'")
    if path is None:
        return open(os.devnull, mode)
    if path == STDIO_PATH:
        if mode == 'r':
            return sys.stdin
        else:
            return sys.stdout
    func = open
    if path.endswith(GZIP_EXT) or path.endswith(BGZIP_EXT):
        func = gzip.open
        mode += 't'
    if path.endswith(BZIP2_EXT):
        func = bz2.open
        mode += 't'
    return func(path, mode)


def count_lines(path):
    with open_file(path) as f:
        count = 0
        for line in f:
            count += 1
    return count


def export_table(table, output_path, header=None):
    with open_file(output_path, 'w') as output_file:
        output_csv = csv.writer(output_file, dialect='excel-tab')
        if header is not None:
            output_csv.writerow(header)
        for row in table:
            output_csv.writerow(row)


# ARGUMENT CHECKS
def flatten_list(initial):
    if not isinstance(initial, (list, tuple)):
        return [initial]
    else:
        final = []
        for i in initial:
            if not isinstance(initial, (list, tuple)):
                final.append(i)
            else:
                final.extend(i)
    return final


def format_errors(errors, prefix=None):
    if prefix is None:
        return errors
    return ['{}: {}'.format(prefix, e) for e in errors]


def check_file_arg(path, mode='r', ext=None, zip_ext=None,
                   none_allowed=False, stdio_allowed=False,
                   prefix=None):
    if mode not in ('r', 'w'):
        raise ValueError("Argument 'mode' must be 'r' or 'w'")
    errors = []
    if path is None:
        if not none_allowed:
            error = 'A path must be specified.'
            errors.append(error)
        return format_errors(errors, prefix)
    if path == STDIO_PATH:
        if not stdio_allowed:
            if mode == 'r':
                error = 'Cannot read from standard input.'
            elif mode == 'w':
                error = 'Cannot write to standard output.'
            errors.append(error)
        return format_errors(errors, prefix)
    if ext is not None:
        ext = flatten_list(ext)
        allowed_extensions = list(ext)
        if zip_ext is not None:
            zip_ext = flatten_list(zip_ext)
            for z in zip_ext:
                allowed_extensions += [e + z for e in ext]

        basename = os.path.basename(path).lower()
        for e in allowed_extensions:
            if basename.endswith(e.lower()):
                break
        else:
            error = 'File {} does not have a good extension.'.format(path)
            errors.append(error)
    try:
        file = open_file(path, mode)
        file.close()
        if mode == 'w':
            os.remove(path)
    except FileNotFoundError:
        error = 'File {} not found.'.format(path)
        errors.append(error)
    except IsADirectoryError:
        error = '{} is a directory.'.format(path)
        errors.append(error)
    except PermissionError:
        permission = 'readable' if mode == 'r' else 'writable'
        error = 'File {} is not {}.'.format(path, permission)
        errors.append(error)
    except (OSError, Exception):
        error = 'An error occurred while opening {}.'.format(path)
        errors.append(error)
    return format_errors(errors, prefix)


def check_dir_arg(path, mode='r', create=False,
                  none_allowed=False, prefix=None):
    if mode not in ('r', 'w'):
        raise ValueError("Argument 'mode' must be 'r' or 'w'")
    errors = []
    if path is None:
        if not none_allowed:
            error = 'A path must be specified.'
            errors.append(error)
        return format_errors(errors, prefix)
    if os.path.exists(path):
        if not os.path.isdir(path):
            errors.append("{}: {} is not a directory.".format(path))
        check = os.R_OK if mode == 'r' else os.W_OK
        permission = 'readable' if mode == 'r' else 'writable'
        if not os.access(path, check):
            errors.append('{}: Directory {} is not '
                          '{}.'.format(path, permission))
    elif mode == 'w' and create:
        try:
            os.mkdir(path)
        except FileNotFoundError:
            errors.append('{}: Directory {} could not be '
                          'created because parent '
                          'directory does not exist.').format(path)
        except PermissionError:
            errors.append('{}: Directory {} could not be '
                          'created because parent '
                          'directory is not writable.'.format(path))
        except Exception:
            errors.append('{}: Directory {} does not exist and could not '
                          'be created.'.format(path))
    return format_errors(errors, prefix)


def check_bin_arg(prog, executable=True, none_allowed=False, prefix=None):
    errors = []
    if prog is None:
        if not none_allowed:
            error = 'A path must be specified.'
            errors.append(error)
        return format_errors(errors, prefix)
    folder, basename = os.path.split(prog)
    if folder:
        if not os.path.exists(prog):
            errors.append('File {} not found.'.format(prog))
        elif not os.path.isfile(prog):
            errors.append('{} is not a file.'.format(prog))
        elif executable and not os.access(prog, os.X_OK):
            errors.append('{} is not executable.'.format(prog))
    else:
        for folder in os.environ['PATH'].split(os.pathsep):
            folder = folder.strip('"')
            path = os.path.join(folder, basename)
            if os.path.exists(path) and os.path.isfile(path):
                if executable and os.access(path, os.X_OK):
                    break
        else:
            errors.append('File {} not found in the '
                          'PATH variable.'.format(path))
    return format_errors(errors, prefix)


def check_num_arg(value, number_type=int, mini=None,
                  maxi=None, none_allowed=False, prefix=None):
    if not isinstance(number_type, type):
        raise ValueError("Argument 'number_type' must be of type 'type'")
    if number_type not in (int, float):
        raise ValueError("Argument 'number_type' must be 'int' or 'float'")
    t = 'int' if number_type is int else 'int or float'
    good_type = int if number_type is int else (int, float)
    if mini is not None and not isinstance(mini, good_type):
        raise TypeError("Argument 'mini' must be of type {}".format(t))
    if maxi is not None and not isinstance(maxi, good_type):
        raise TypeError("Argument 'maxi' must be of type {}".format(t))
    if mini is not None and maxi is not None and mini > maxi:
        raise ValueError("Argument 'mini' must <= argument 'maxi'")
    errors = []
    if value is None:
        if not none_allowed:
            error = 'A value must be specified.'
            errors.append(error)
        return format_errors(errors, prefix)
    if not isinstance(value, number_type):
        error = "Value '{}' is not of type {}.".format(value, t)
        errors.append(error)
    else:
        if mini is not None and maxi is not None:
            if mini == maxi and value != mini:
                error = "Value must be equal to {!s}.".format(mini)
                errors.append(error)
            elif not mini <= value <= maxi:
                error = ("Value must be chosen between {!s} "
                         "and {!s}.".format(mini, maxi))
                errors.append(error)
        elif mini is not None and value < mini:
            error = "Value must be >= {!s}.".format(mini)
            errors.append(error)
        elif maxi is not None and value > maxi:
            error = "Value must be <= {!s}.".format(maxi)
            errors.append(error)
    return format_errors(errors, prefix)


def check_time_arg(value, time_format='%H:%M', none_allowed=False,
                   prefix=None):
    errors = []
    if value is None:
        if not none_allowed:
            error = 'A value must be specified.'
            errors.append(error)
        return format_errors(errors, prefix)
    try:
        datetime.datetime.strptime(value, time_format)
    except ValueError:
        error = "Wrong time format for value '{}'.".format(value)
        errors.append(error)
    return format_errors(errors, prefix)


def check_seq_file_format(path, specified_format,
                          default_format=SEQFILE_FORMAT,
                          file_indexed=FILE_INDEXED):
    if path == STDIO_PATH:
        if specified_format is not None:
            return specified_format
        return default_format
    basename = os.path.basename(path)
    for ext in ZIP_EXTS:
        if not basename.endswith(ext):
            continue
        unzipped_basename, zip_ext = basename[:-1*len(ext)], ext
        break
    else:
        unzipped_basename, zip_ext = basename, None

    bgzip_error = (
        'If an index file is provided, the sequence file must '
        'be uncompressed, or compressed in BGZ format.')
    wrong_ext_error = (
        'File extension does not correspond to the specified format.')
    unknown_ext_error = (
        'File extension does not correspond to any known format.')

    for f, exts in zip([FASTA_FORMAT, GENBANK_FORMAT],
                       [FASTA_EXTS, GENBANK_EXTS]):
        if specified_format in [None, f]:
            for e in exts:
                if unzipped_basename.endswith(e):
                    if not file_indexed:
                        return f
                    if zip_ext is None:
                        return f
                    if zip_ext == BGZIP_EXT:
                        return f
                    raise ValueError(bgzip_error)
            if specified_format == f:
                raise ValueError(wrong_ext_error)
    raise ValueError(unknown_ext_error)


# SUBPROCESSING
def run_child_process(cmd, stdout_path=None, stderr_path=None,
                      timeout=None, remove_log=False,
                      exit_if_error=False, print_msg=False,
                      program=None):
    cmd = [str(x) for x in cmd]
    if program is None:
        program = os.path.basename(cmd[0])
    if stdout_path is None:
        stdout_path = program.replace(' ', '_') + '.log'
    if stderr_path is None:
        stderr_path = program.replace(' ', '_') + '.err'
    merged_outputs = True
    if stdout_path != stderr_path:
        merged_outputs = False
    if print_msg:
        print_stderr('Launching command...', program=program)
    try:
        stdout = open(stdout_path, 'w')
        if not merged_outputs:
            stderr = open(stderr_path, 'w')
        else:
            stderr = stdout
        p = subprocess.run(cmd, stdout=stdout, stderr=stderr,
                           timeout=timeout, check=True)
        stdout.close()
        if not merged_outputs:
            stderr.close()
    except subprocess.CalledProcessError as exc:
        if print_msg:
            cmd_str = ' '.join(cmd)
            msg = ('Process return a non-zero exit code ({})\n'
                   'See file: {}\n'
                   'Command: {}').format(exc.returncode, stderr_path, cmd_str)
            print_stderr(msg, program=program, prefix='ERROR')
        if exit_if_error:
            sys.exit(1)
    except subprocess.TimeoutExpired as exc:
        if print_msg:
            cmd_str = ' '.join(cmd)
            msg = ('Process timed out ({})\n'
                   'See file: {}\n'
                   'Command: {}').format(exc.timeout, stderr_path, cmd_str)
            print_stderr(msg, program=program, prefix='ERROR')
        if exit_if_error:
            sys.exit(1)
    else:
        if remove_log:
            os.remove(stdout_path)
            if stderr_path != stdout_path:
                os.remove(stderr_path)
        if print_msg:
            print_stderr('Done.', program=program, end='\n\n')
    return p


# SEQUENCE I/O
def import_sequences(input_path, input_format=FASTA_FORMAT):
    return Bio.SeqIO.to_dict(Bio.SeqIO.parse(input_path, input_format))


def export_sequences(records, output_path, output_format=FASTA_FORMAT):
    if isinstance(records, dict):
        records = [records[seqid] for seqid in sorted(records.keys())]
    Bio.SeqIO.write(records, output_path, output_format)


def write_records(records, output_file, output_format=FASTA_FORMAT):
    if isinstance(records, Bio.SeqRecord.SeqRecord):
        records = [records]
    override_alphabet = False
    if output_format == GENBANK_FORMAT:
        output_alphabet = find_records_alphabet(*records)
        if not isinstance(output_alphabet, (NT_ALPHA, AA_ALPHA)):
            override_alphabet = True
            output_alphabet = guess_records_alphabet(*records)
            if output_alphabet is None:
                raise ValueError('Sequence alphabet could not be determined.')
    for record in records:
        if override_alphabet:
            record.seq.alphabet = output_alphabet
        Bio.SeqIO.write(record, output_file, output_format)


def find_records_alphabet(*records, sample_size=100):
    alphabet = None
    nb_rec = 0
    for record in records:
        nb_rec += 1
        if nb_rec > sample_size:
            return alphabet
        rec_alphabet = record.seq.alphabet
        if rec_alphabet is None:
            return None
        if alphabet is not None and rec_alphabet != alphabet:
            return None
        rec_alphabet = alphabet
    return alphabet


def guess_records_alphabet(*records, sample_size=10000):
    sample = ''
    for record in records:
        sample += str(record.seq)
        if len(sample) < sample_size:
            break

    characters = set(list(sample.upper())) - {'*', '-', '.'}

    if characters < set(DNA_ALPHA_EXT.letters):
        if characters == set(DNA_ALPHA.letters):
            return DNA_ALPHA()
        return DNA_ALPHA_EXT()

    if characters < set(PROT_ALPHA_EXT.letters):
        if characters == set(PROT_ALPHA.letters):
            return PROT_ALPHA()
        return PROT_ALPHA_EXT()

    if characters < set(RNA_ALPHA_EXT.letters):
        if characters == set(RNA_ALPHA.letters):
            return RNA_ALPHA()
        return RNA_ALPHA_EXT()

    return None


def decode_gff_attr(attr_str, sep='='):
    """Converts GFF/GTF attribute string to a dictionary"""
    attr = {}
    pattern = re.compile('''((?:[^;"']|"[^"]*"|'[^']*')+)''')
    for key_value in pattern.findall(attr_str):
        key_value = key_value.strip()
        try:
            pos = key_value.index(sep)
            key, value = key_value[:pos], key_value[pos+1:]
        except:
            print(attr_str, key_value)
        attr[key] = value.strip('"')
    return attr


def decode_gtf_attr(attr_str):
    return decode_gff_attr(attr_str, sep=' ')


# SEQUENCE OPERATIONS
def find_consensus_nt(nts):
    ambiguous_nts = {
        ('A',): 'A',
        ('C',): 'C',
        ('G',): 'G',
        ('T',): 'T',
        ('A', 'C'): 'M',
        ('A', 'G'): 'R',
        ('A', 'T'): 'W',
        ('C', 'G'): 'S',
        ('C', 'T'): 'Y',
        ('G', 'T'): 'K',
        ('A', 'C', 'G'): 'V',
        ('A', 'C', 'T'): 'H',
        ('A', 'G', 'T'): 'D',
        ('C', 'G', 'T'): 'B'}
    nts = tuple(sorted(set(nts)))
    try:
        return ambiguous_nts[nts]
    except KeyError:
        return 'N'


# PARTITIONING
def partition_list(values, nb_parts=1):
    nb_values = len(values)
    if nb_values <= nb_parts:
        return [[i] for i in values]
    parts = []
    part_length = nb_values // nb_parts
    modulo = nb_values % nb_parts
    start = 0
    for n in range(nb_parts):
        end = start + part_length
        if n < modulo:
            end += 1
        parts.append([values[i] for i in range(start, end)])
        start = end
    return parts


def merge_dictionaries(dicts, merge_func=None):
    if dicts is None or not len(dicts):
        return dicts
    if len(dicts) == 1:
        return dicts[0]
    merged_dict = {}
    for d in dicts:
        for key in d:
            try:
                merged_dict[key].append(d[key])
            except KeyError:
                merged_dict[key] = [d[key]]
    if merge_func is None:
        return merged_dict
    for key in merged_dict:
        merged_dict[key] = merge_func(merged_dict[key])
    return merged_dict


def updict_append_to_list(d, key, item):
    """Update a dict of lists by appending a new item.

    If 'key' is in 'd', append 'item' to 'd[key]'. If 'key' is not in
    'd', create a entry 'd[key]' containing '[item]'.

    Arguments:
        d (dict): An existing dictionary.
        key: The dictionary key.
        item: The new item.

    Raises:
        TypeError: If 'd' is not a dictionary or 'key' is not hashable.
    """
    if not isinstance(d, dict):
        raise TypeError("Argument 'd' must be a dictionary.")
    if not isinstance(key, collections.Hashable):
        raise TypeError("Argument 'key' must be hashable.")
    try:
        d[key].append(item)
    except KeyError:
        d[key] = [item]


def updict_extend_list(d, key, new_list):
    """Update a dict of lists by extending an existing list.

    If 'key' is in 'd', append 'new_list' to 'd[key]'. If 'key' is not
        in 'd', create a new entry in 'd[key]' containing 'new_list'.

    Arguments:
        d (dict): An existing dictionary.
        key (hashable): The dictionary key.
        item: The new item.

    Raises:
        TypeError: If 'd' is not a dictionary, 'key' is not hashable, or
            'new_list' is not a list, tuple, set, or frozenset.
    """
    if not isinstance(d, dict):
        raise TypeError("Argument 'd' must be a dictionary.")
    if not isinstance(key, collections.Hashable):
        raise TypeError("Argument 'key' must be hashable.")
    if not isinstance(new_list, (list, tuple, set, frozenset)):
        raise TypeError("Argument 'new_list' must be a list, tuple, "
                        "set, or frozenset.")
    try:
        d[key].extend(new_list)
    except KeyError:
        d[key] = list(new_list)


def updict_add_to_set(d, key, item):
    """Update a dict of sets by adding a new item to an existing set.

    If 'key' is in 'd', append 'new_list' to 'd[key]'. If 'key' is not
        in 'd', create a new list in 'd[key]' containing '[new_list]'.

    Arguments:
        d (dict): An existing dictionary.
        key (hashable): The dictionary key.
        new_list (list, tuple, set, frozenset): The new iterable used to
            extend the existing list.

    Raises:
        TypeError: If 'd' is not a dictionary or 'key' is not hashable.
    """
    if not isinstance(d, dict):
        raise TypeError("Argument 'd' must be a dictionary.")
    if not isinstance(key, collections.Hashable):
        raise TypeError("Argument 'key' must be hashable.")
    try:
        d[key] |= {item}
    except KeyError:
        d[key] = {item}


def updict_set_union(d, key, new_set):
    """Update a dict of sets by extending an existing set.

    If 'key' is in 'd', append 'new_list' to 'd[key]'. If 'key' is not
        in 'd', create a new entry in 'd[key]' containing 'new_list'.

    Arguments:
        d (dict): An existing dictionary.
        key (hashable): The dictionary key.
        new_set (set, frozenset): The new set used to extend the
            existing set.

    Raises:
        TypeError: If 'd' is not a dictionary, 'key' is not hashable, or
            'new_set' is not a set nor a frozenset.
    """
    if not isinstance(d, dict):
        raise TypeError("Argument 'd' must be a dictionary.")
    if not isinstance(key, collections.Hashable):
        raise TypeError("Argument 'key' must be hashable.")
    if not isinstance(new_set, (set, frozenset)):
        raise TypeError("Argument 'new_set' must be a set or frozenset.")
    try:
        d[key] |= new_set
    except KeyError:
        d[key] = set(new_set)
