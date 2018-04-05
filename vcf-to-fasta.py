# Requires pyVCF
# The order of bases in each sequence may vary depending on python version.
# Tested on python 2.7 and 3.6
from __future__ import print_function
import sys
from collections import OrderedDict
from re import split
import random
import argparse
import vcf  # TODO make it work without pyVCF

INCLUDE_REFERENCE = False
USE_IUPAC_AMBIGUITY_CODE = False
QUALITY_THRESHOLD = 20
TREAT_LOWQUALITY = ['as_genotype', 'stochastic', 'ambiguous', 'delete'][0]
# TREAT_MISSING = ['as_reference', 'as_N', 'delete'][0]

IUPAC_CODE = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
              'AG': 'R', 'CT': 'Y', 'GT': 'K', 'AC': 'M', 'CG': 'S', 'AT': 'W',
              'CGT': 'B', 'AGT': 'D', 'ACT': 'H', 'ACG': 'V'}

def cat_gt_bases(gt_bases):
    return ''.join(split('[/|]', gt_bases))

def iupac(bases):
    ambigcode = IUPAC_CODE.get(''.join(sorted(set(bases))), default='N')
    return ambigcode

def stochastic_selection(rec, bases):
    p = 10**-(rec.QUAL/10)
    if random.random() > p:
        return bases
    else:
        return rec.REF


class VcfConflictError(Exception):
    pass

def extract_snp_from_vcf(*vcffiles):
    snps_data = OrderedDict()  # sample => dict of records

    contigs_foolproof = None
    for vcffile in vcffiles:
        try:
            vcf_reader = vcf.Reader(filename=vcffile)  # takes .vcf or .vcf.gz file
        except UnicodeDecodeError:
            print('ERROR: Invalid input VCF file(s)', file=sys.stderr)
            sys.exit(1)

        if not contigs_foolproof:
            contigs_foolproof = vcf_reader.contigs
        elif vcf_reader.contigs != contigs_foolproof:
            print('ERROR: Different contigs information between multiple vcf files', file=sys.stderr)
            raise VcfConflictError

        if 'unknown' in vcf_reader.samples:
            print('WARNING: Sample name not designated in vcf file:', vcffile, file=sys.stderr)

        for rec in vcf_reader:
            if not rec.is_snp:  # if the record is indel or complex
                continue

            if rec.QUAL < QUALITY_THRESHOLD:
                print('WARNING: Low quality record at {0} from {1}'.format(rec, vcffile), 
                      file=sys.stderr)
                is_rec_low_qual = True
            else:
                is_rec_low_qual = False

            try:
                # if there is a conflict between references from multiple records
                if rec.REF != snps_data.setdefault('reference', {})[rec.POS]:  # (initialize dict at the same time)
                    print('ERROR: Multiple reference bases in the same position {}'.format(rec.POS),
                          file=sys.stderr)
                    raise VcfConflictError
            except KeyError:
                pass
            finally:
                snps_data['reference'][rec.POS] = rec.REF  # set reference base at the position

            for call in rec.samples:
                if call.sample == 'unknown':
                    sample = vcf_reader.filename  # ! only a makeshift; may cause sample name confliction
                else:
                    sample = call.sample

                if not call.called:  # missing call as 'N'
                    snps_data.setdefault(sample, {})[rec.POS] = 'N'  # TODO
                    continue

                bases = cat_gt_bases(call.gt_bases)  # genotype bases for the call
                if handle_lowquality and is_rec_low_qual:
                    bases = handle_lowquality(rec, bases)

                if not call.is_het:  # if homozygous
                    base = bases[0]
                else:
                    base = handle_heterozygous(bases)
                snps_data.setdefault(sample, {})[rec.POS] = base

    return snps_data


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('vcffiles', metavar='<VCF file>', nargs='+')
    parser.add_argument('--include-reference', dest='incl_ref', action='store_true',
                        help='Include reference sequence in the final output')
    parser.add_argument('--iupac', action='store_true',
                        help='For heterozygous SNPs, use IUPAC ambiguity code instead of random sampling')
    parser.add_argument('--quality', metavar='X', dest='qualchk', type=float, default=20,
                        help='Warn when a record has QUAL score less than X (0 to disable warning. Default 20)')
    parser.add_argument('--treat-low-qual', dest='treatlowq',
                        choices=['as_genotype', 'as_reference', 'as_ambiguous', 'stochastic'], default='as_genotype',
                        help='##########')  # TODO description
    # parser.add_argument('--treat-missing', dest='treatmiss',
    #                     choices=['as_reference', 'as_N', 'delete'], default='as_reference',
    #                     help='##########')
    args = parser.parse_args()

    vcffiles = args.vcffiles
    INCLUDE_REFERENCE = args.incl_ref
    USE_IUPAC_AMBIGUITY_CODE = args.iupac
    QUALITY_THRESHOLD = args.qualchk
    TREAT_LOWQUALITY = args.treatlowq
    # TREAT_MISSING = args.treatmiss

    if USE_IUPAC_AMBIGUITY_CODE:
        handle_heterozygous = iupac
    else:
        handle_heterozygous = random.choice  # random haplotyping

    handle_lowquality = {
        'as_reference': lambda rec, _: rec.REF,
        'as_ambiguous': lambda _, bases: iupac(bases),
        'stochastic': stochastic_selection
    }.get(TREAT_LOWQUALITY, None)

    try:
        snps_data = extract_snp_from_vcf(*vcffiles)
    except ValueError as e:
        print('ERROR: VCF parse error - maybe wrong VCF format?', file=sys.stderr)
        sys.exit(1)
    except VcfConflictError:
        sys.exit(1)

    columns = []
    for position in snps_data['reference']:
        column = [pos_alt.get(position, snps_data['reference'][position]) for pos_alt in snps_data.values()]
        first_elem = column[not INCLUDE_REFERENCE]
        if any(x != first_elem for x in column):
            columns.append(column)

    en = enumerate(snps_data)
    if not INCLUDE_REFERENCE:
        next(en)
    for i, sample in en:
        print('>{}'.format(sample))
        print(''.join(column[i] for column in columns))
