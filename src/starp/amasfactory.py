from .parsers import TwoAlleles
from .models import Sequence, Snp, AmasPrimer
from .utils import hamming

# These substitution matrices were created to deal with the
# combinatorial explosion of possible primer endings given by
# docs/starp/STARP F primer design_clarified.docx
#
# Exempting the final nucleotide, 'G' and 'C' nucleotides were changed
# to an 'S' while 'A' and 'T' nucleotides were changed to a 'W'.
# 'P' represents a [C/G] OR [A/T] snps while 'N' represents all
# other substitution snps.

# S -> G or C
# W -> A or T
# P -> G/C or A/T SNP
# N -> A/G, A/C, T/G, or T/C SNP

sub_index_one_snp = {
    frozenset(('C', 'G')) : {
        'SSSC' : -3, 'SSSG' : -2, 'SSWC' : -4, 'SSWG' : -3,
        'SWSC' : -4, 'SWSG' : -2, 'SWWC' : -3, 'SWWG' : -2,
        'WSSC' : -3, 'WSSG' : -2, 'WSWC' : -4, 'WSWG' : -2,
        'WWSC' : -4, 'WWSG' : -3, 'WWWC' : -3, 'WWWG' : -2},
    frozenset(('C', 'T')) : {
        'SSSC' : -2, 'SSST' : -3, 'SSWC' : -4, 'SSWT' : -3,
        'SWSC' : -2, 'SWST' : -4, 'SWWC' : -2, 'SWWT' : -3,
        'WSSC' : -2, 'WSST' : -3, 'WSWC' : -2, 'WSWT' : -4,
        'WWSC' : -4, 'WWST' : -3, 'WWWC' : -2, 'WWWT' : -3},
    frozenset(('C', 'A')) : {
        'SSSC' : -2, 'SSSA' : -3, 'SSWC' : -3, 'SSWA' : -4,
        'SWSC' : -2, 'SWSA' : -4, 'SWWC' : -2, 'SWWA' : -3,
        'WSSC' : -2, 'WSSA' : -3, 'WSWC' : -2, 'WSWA' : -4,
        'WWSC' : -3, 'WWSA' : -4, 'WWWC' : -2, 'WWWA' : -3},
    frozenset(('G', 'T')) : {
        'SSSG' : -2, 'SSST' : -3, 'SSWG' : -4, 'SSWT' : -3,
        'SWSG' : -2, 'SWST' : -4, 'SWWG' : -2, 'SWWT' : -3,
        'WSSG' : -2, 'WSST' : -3, 'WSWG' : -2, 'WSWT' : -4,
        'WWSG' : -4, 'WWST' : -3, 'WWWG' : -2, 'WWWT' : -3},
    frozenset(('G', 'A')) : {
        'SSSG' : -2, 'SSSA' : -3, 'SSWG' : -3, 'SSWA' : -4,
        'SWSG' : -2, 'SWSA' : -4, 'SWWG' : -2, 'SWWA' : -3,
        'WSSG' : -2, 'WSSA' : -3, 'WSWG' : -2, 'WSWA' : -4,
        'WWSG' : -3, 'WWSA' : -4, 'WWWG' : -2, 'WWWA' : -3},
    frozenset(('A', 'T')) : {
        'SSST' : -3, 'SSSA' : -4, 'SSWT' : -3, 'SSWA' : -4,
        'SWST' : -2, 'SWSA' : -4, 'SWWT' : -2, 'SWWA' : -3,
        'WSST' : -2, 'WSSA' : -3, 'WSWT' : -2, 'WSWA' : -4,
        'WWST' : -3, 'WWSA' : -4, 'WWWT' : -3, 'WWWA' : -4}}

sub_index_two_snps = {
    frozenset(('C', 'G')) : {

        # 1.1 A
        'WWWPC': {'S': -4, 'W': -4}, 'WWWPG': {'S': -3, 'W': -3}, 
        'WWSPC': {'S': -5, 'W': -5}, 'WWSPG': {'S': -4, 'W': -4},
        'WSWPC': {'S': -5, 'W': -5}, 'WSWPG': {'S': -3, 'W': -3},
        'SWWPC': {'S': -4, 'W': -4}, 'SWWPG': {'S': -3, 'W': -3},
        'WSSPC': {'S': -4, 'W': -4}, 'WSSPG': {'S': -3, 'W': -3},
        'SWSPC': {'S': -5, 'W': -5}, 'SWSPG': {'S': -3, 'W': -3},
        'SSWPC': {'S': -5, 'W': -5}, 'SSWPG': {'S': -4, 'W': -4},
        'SSSPC': {'S': -4, 'W': -4}, 'SSSPG': {'S': -3, 'W': -3},

        # 1.1 B
        'WWWNC': {'S': -4, 'W': -4}, 'WWWNG': {'S': -3, 'W': -3}, 
        'WWSNC': {'S': -3, 'W': -4}, 'WWSNG': {'S': -3, 'W': -4},
        'WSWNC': {'S': -4, 'W': -3}, 'WSWNG': {'S': -4, 'W': -3},
        'SWWNC': {'S': -4, 'W': -4}, 'SWWNG': {'S': -3, 'W': -3},
        'WSSNC': {'S': -4, 'W': -4}, 'WSSNG': {'S': -3, 'W': -3},
        'SWSNC': {'S': -3, 'W': -4}, 'SWSNG': {'S': -3, 'W': -4},
        'SSWNC': {'S': -4, 'W': -3}, 'SSWNG': {'S': -4, 'W': -3},
        'SSSNC': {'S': -4, 'W': -4}, 'SSSNG': {'S': -3, 'W': -3},

        # 1.2 A
        'WWPWC': {'S': -4, 'W': -4}, 'WWPWG': {'S': -2, 'W': -2}, 
        'WWPSC': {'S': -5, 'W': -5}, 'WWPSG': {'S': -4, 'W': -4},
        'WSPWC': {'S': -5, 'W': -5}, 'WSPWG': {'S': -2, 'W': -2},
        'SWPWC': {'S': -4, 'W': -4}, 'SWPWG': {'S': -2, 'W': -2},
        'WSPSC': {'S': -4, 'W': -4}, 'WSPSG': {'S': -2, 'W': -2},
        'SWPSC': {'S': -5, 'W': -5}, 'SWPSG': {'S': -2, 'W': -2},
        'SSPWC': {'S': -5, 'W': -5}, 'SSPWG': {'S': -4, 'W': -4},
        'SSPSC': {'S': -4, 'W': -4}, 'SSPSG': {'S': -2, 'W': -2},

        # 1.2 B
        'WWNWC': {'S': -4, 'W': -4}, 'WWNWG': {'S': -2, 'W': -2}, 
        'WWNSC': {'S': -2, 'W': -4}, 'WWNSG': {'S': -2, 'W': -4},
        'WSNWC': {'S': -4, 'W': -2}, 'WSNWG': {'S': -4, 'W': -2},
        'SWNWC': {'S': -4, 'W': -4}, 'SWNWG': {'S': -2, 'W': -2},
        'WSNSC': {'S': -4, 'W': -4}, 'WSNSG': {'S': -2, 'W': -2},
        'SWNSC': {'S': -2, 'W': -4}, 'SWNSG': {'S': -2, 'W': -4},
        'SSNWC': {'S': -4, 'W': -2}, 'SSNWG': {'S': -4, 'W': -2},
        'SSNSC': {'S': -4, 'W': -4}, 'SSNSG': {'S': -2, 'W': -2},

        # 1.3 A
        'WPWWC': {'S': -3, 'W': -3}, 'WPWWG': {'S': -2, 'W': -2}, 
        'WPWSC': {'S': -5, 'W': -5}, 'WPWSG': {'S': -3, 'W': -3},
        'WPSWC': {'S': -5, 'W': -5}, 'WPSWG': {'S': -2, 'W': -2},
        'SPWWC': {'S': -3, 'W': -3}, 'SPWWG': {'S': -2, 'W': -2},
        'WPSSC': {'S': -3, 'W': -3}, 'WPSSG': {'S': -2, 'W': -2},
        'SPWSC': {'S': -5, 'W': -5}, 'SPWSG': {'S': -2, 'W': -2},
        'SPSWC': {'S': -5, 'W': -5}, 'SPSWG': {'S': -3, 'W': -3},
        'SPSSC': {'S': -3, 'W': -3}, 'SPSSG': {'S': -2, 'W': -2},

        # 1.3 B
        'WNWWC': {'S': -3, 'W': -3}, 'WNWWG': {'S': -2, 'W': -2},
        'WNWSC': {'S': -2, 'W': -3}, 'WNWSG': {'S': -2, 'W': -3},
        'WNSWC': {'S': -3, 'W': -2}, 'WNSWG': {'S': -3, 'W': -2},
        'SNWWC': {'S': -3, 'W': -3}, 'SNWWG': {'S': -2, 'W': -2},
        'WNSSC': {'S': -3, 'W': -3}, 'WNSSG': {'S': -2, 'W': -2},
        'SNWSC': {'S': -2, 'W': -3}, 'SNWSG': {'S': -2, 'W': -3},
        'SNSWC': {'S': -3, 'W': -2}, 'SNSWG': {'S': -3, 'W': -2},
        'SNSSC': {'S': -3, 'W': -3}, 'SNSSG': {'S': -2, 'W': -2},

        # 1.4 A
        'PWWWC': {'S': -3, 'W': -3}, 'PWWWG': {'S': -2, 'W': -2}, 
        'PWWSC': {'S': -4, 'W': -4}, 'PWWSG': {'S': -3, 'W': -3},
        'PWSWC': {'S': -4, 'W': -4}, 'PWSWG': {'S': -2, 'W': -2},
        'PSWWC': {'S': -3, 'W': -3}, 'PSWWG': {'S': -2, 'W': -2},
        'PWSSC': {'S': -3, 'W': -3}, 'PWSSG': {'S': -2, 'W': -2},
        'PSWSC': {'S': -4, 'W': -4}, 'PSWSG': {'S': -2, 'W': -2},
        'PSSWC': {'S': -4, 'W': -4}, 'PSSWG': {'S': -3, 'W': -3},
        'PSSSC': {'S': -3, 'W': -3}, 'PSSSG': {'S': -2, 'W': -2},

        # 1.4 B
        'NWWWC': {'S': -3, 'W': -3}, 'NWWWG': {'S': -2, 'W': -2},
        'NWWSC': {'S': -2, 'W': -3}, 'NWWSG': {'S': -2, 'W': -3},
        'NWSWC': {'S': -3, 'W': -2}, 'NWSWG': {'S': -3, 'W': -2},
        'NSWWC': {'S': -3, 'W': -3}, 'NSWWG': {'S': -2, 'W': -2},
        'NWSSC': {'S': -3, 'W': -3}, 'NWSSG': {'S': -2, 'W': -2},
        'NSWSC': {'S': -2, 'W': -3}, 'NSWSG': {'S': -2, 'W': -3},
        'NSSWC': {'S': -3, 'W': -2}, 'NSSWG': {'S': -3, 'W': -2},
        'NSSSC': {'S': -3, 'W': -3}, 'NSSSG': {'S': -2, 'W': -2},

    },

    frozenset(('C', 'T')) : {

        # 2.1 A
        'WWWPC': {'S': -4, 'W': -4}, 'WWWPT': {'S': -3, 'W': -3}, 
        'WWSPC': {'S': -3, 'W': -3}, 'WWSPT': {'S': -4, 'W': -4},
        'WSWPC': {'S': -4, 'W': -4}, 'WSWPT': {'S': -3, 'W': -3},
        'SWWPC': {'S': -4, 'W': -4}, 'SWWPT': {'S': -3, 'W': -3},
        'WSSPC': {'S': -4, 'W': -4}, 'WSSPT': {'S': -3, 'W': -3},
        'SWSPC': {'S': -3, 'W': -3}, 'SWSPT': {'S': -3, 'W': -3},
        'SSWPC': {'S': -4, 'W': -4}, 'SSWPT': {'S': -4, 'W': -4},
        'SSSPC': {'S': -4, 'W': -4}, 'SSSPT': {'S': -3, 'W': -3},

        # 2.1 B
        'WWWNC': {'S': -4, 'W': -4}, 'WWWNT': {'S': -3, 'W': -3}, 
        'WWSNC': {'S': -3, 'W': -5}, 'WWSNT': {'S': -3, 'W': -4},
        'WSWNC': {'S': -4, 'W': -4}, 'WSWNT': {'S': -3, 'W': -3},
        'SWWNC': {'S': -5, 'W': -4}, 'SWWNT': {'S': -3, 'W': -3},
        'WSSNC': {'S': -3, 'W': -5}, 'WSSNT': {'S': -3, 'W': -5},
        'SWSNC': {'S': -3, 'W': -5}, 'SWSNT': {'S': -4, 'W': -4},
        'SSWNC': {'S': -4, 'W': -4}, 'SSWNT': {'S': -3, 'W': -3},
        'SSSNC': {'S': -4, 'W': -5}, 'SSSNT': {'S': -4, 'W': -4},

        # 2.2 A
        'WWPWC': {'S': -2, 'W': -2}, 'WWPWT': {'S': -4, 'W': -4}, 
        'WWPSC': {'S': -2, 'W': -2}, 'WWPST': {'S': -4, 'W': -4},
        'WSPWC': {'S': -2, 'W': -2}, 'WSPWT': {'S': -5, 'W': -5},
        'SWPWC': {'S': -2, 'W': -2}, 'SWPWT': {'S': -4, 'W': -4},
        'WSPSC': {'S': -2, 'W': -2}, 'WSPST': {'S': -4, 'W': -4},
        'SWPSC': {'S': -2, 'W': -2}, 'SWPST': {'S': -4, 'W': -4},
        'SSPWC': {'S': -5, 'W': -5}, 'SSPWT': {'S': -4, 'W': -4},
        'SSPSC': {'S': -2, 'W': -2}, 'SSPST': {'S': -4, 'W': -4},

        # 2.2 B
        'WWNWC': {'S': -2, 'W': -2}, 'WWNWT': {'S': -4, 'W': -4}, 
        'WWNSC': {'S': -2, 'W': -5}, 'WWNST': {'S': -4, 'W': -4},
        'WSNWC': {'S': -4, 'W': -2}, 'WSNWT': {'S': -5, 'W': -5},
        'SWNWC': {'S': -5, 'W': -2}, 'SWNWT': {'S': -4, 'W': -4},
        'WSNSC': {'S': -2, 'W': -2}, 'WSNST': {'S': -4, 'W': -5},
        'SWNSC': {'S': -2, 'W': -2}, 'SWNST': {'S': -5, 'W': -4},
        'SSNWC': {'S': -5, 'W': -5}, 'SSNWT': {'S': -4, 'W': -4},
        'SSNSC': {'S': -2, 'W': -2}, 'SSNST': {'S': -4, 'W': -4},

        # 2.3 A
        'WPWWC': {'S': -2, 'W': -2}, 'WPWWT': {'S': -3, 'W': -3}, 
        'WPWSC': {'S': -2, 'W': -2}, 'WPWST': {'S': -3, 'W': -3},
        'WPSWC': {'S': -2, 'W': -2}, 'WPSWT': {'S': -5, 'W': -5},
        'SPWWC': {'S': -2, 'W': -2}, 'SPWWT': {'S': -3, 'W': -3},
        'WPSSC': {'S': -2, 'W': -2}, 'WPSST': {'S': -3, 'W': -3},
        'SPWSC': {'S': -2, 'W': -2}, 'SPWST': {'S': -3, 'W': -3},
        'SPSWC': {'S': -5, 'W': -5}, 'SPSWT': {'S': -3, 'W': -3},
        'SPSSC': {'S': -2, 'W': -2}, 'SPSST': {'S': -3, 'W': -3},

        # 2.3 B
        'WNWWC': {'S': -2, 'W': -2}, 'WNWWT': {'S': -3, 'W': -3},
        'WNWSC': {'S': -2, 'W': -5}, 'WNWST': {'S': -3, 'W': -3},
        'WNSWC': {'S': -3, 'W': -2}, 'WNSWT': {'S': -5, 'W': -5},
        'SNWWC': {'S': -5, 'W': -2}, 'SNWWT': {'S': -3, 'W': -3},
        'WNSSC': {'S': -2, 'W': -2}, 'WNSST': {'S': -3, 'W': -5},
        'SNWSC': {'S': -2, 'W': -2}, 'SNWST': {'S': -5, 'W': -3},
        'SNSWC': {'S': -5, 'W': -5}, 'SNSWT': {'S': -3, 'W': -3},
        'SNSSC': {'S': -2, 'W': -2}, 'SNSST': {'S': -3, 'W': -3},

        # 2.4 A
        'PWWWC': {'S': -2, 'W': -2}, 'PWWWT': {'S': -3, 'W': -3},
        'PWWSC': {'S': -2, 'W': -2}, 'PWWST': {'S': -3, 'W': -3},
        'PWSWC': {'S': -2, 'W': -2}, 'PWSWT': {'S': -4, 'W': -4},
        'PSWWC': {'S': -2, 'W': -2}, 'PSWWT': {'S': -3, 'W': -3},
        'PWSSC': {'S': -2, 'W': -2}, 'PWSST': {'S': -3, 'W': -3},
        'PSWSC': {'S': -2, 'W': -2}, 'PSWST': {'S': -3, 'W': -3},
        'PSSWC': {'S': -3, 'W': -3}, 'PSSWT': {'S': -3, 'W': -3},
        'PSSSC': {'S': -2, 'W': -2}, 'PSSST': {'S': -3, 'W': -3},

        # 2.4 B
        'NWWWC': {'S': -2, 'W': -2}, 'NWWWT': {'S': -3, 'W': -3},
        'NWWSC': {'S': -2, 'W': -4}, 'NWWST': {'S': -3, 'W': -3},
        'NWSWC': {'S': -3, 'W': -2}, 'NWSWT': {'S': -4, 'W': -4},
        'NSWWC': {'S': -4, 'W': -2}, 'NSWWT': {'S': -3, 'W': -3},
        'NWSSC': {'S': -2, 'W': -2}, 'NWSST': {'S': -3, 'W': -4},
        'NSWSC': {'S': -2, 'W': -2}, 'NSWST': {'S': -4, 'W': -3},
        'NSSWC': {'S': -4, 'W': -4}, 'NSSWT': {'S': -3, 'W': -3},
        'NSSSC': {'S': -2, 'W': -2}, 'NSSST': {'S': -3, 'W': -3},

    },
        
    frozenset(('C', 'A')) : {

        # 3.1 A
        'WWWPC': {'S': -3, 'W': -3}, 'WWWPA': {'S': -4, 'W': -4}, 
        'WWSPC': {'S': -3, 'W': -3}, 'WWSPA': {'S': -4, 'W': -4},
        'WSWPC': {'S': -4, 'W': -4}, 'WSWPA': {'S': -3, 'W': -3},
        'SWWPC': {'S': -3, 'W': -3}, 'SWWPA': {'S': -4, 'W': -4},
        'WSSPC': {'S': -3, 'W': -3}, 'WSSPA': {'S': -4, 'W': -4},
        'SWSPC': {'S': -3, 'W': -3}, 'SWSPA': {'S': -4, 'W': -4},
        'SSWPC': {'S': -4, 'W': -4}, 'SSWPA': {'S': -3, 'W': -3},
        'SSSPC': {'S': -3, 'W': -3}, 'SSSPA': {'S': -4, 'W': -4},

        # 3.1 B
        'WWWNC': {'S': -3, 'W': -3}, 'WWWNA': {'S': -4, 'W': -4}, 
        'WWSNC': {'S': -3, 'W': -4}, 'WWSNA': {'S': -5, 'W': -4},
        'WSWNC': {'S': -4, 'W': -3}, 'WSWNA': {'S': -5, 'W': -3},
        'SWWNC': {'S': -5, 'W': -3}, 'SWWNA': {'S': -4, 'W': -3},
        'WSSNC': {'S': -3, 'W': -3}, 'WSSNA': {'S': -4, 'W': -5},
        'SWSNC': {'S': -3, 'W': -3}, 'SWSNA': {'S': -5, 'W': -4},
        'SSWNC': {'S': -4, 'W': -4}, 'SSWNA': {'S': -5, 'W': -3},
        'SSSNC': {'S': -3, 'W': -3}, 'SSSNA': {'S': -4, 'W': -4},

        # 3.2 A
        'WWPWC': {'S': -2, 'W': -2}, 'WWPWA': {'S': -4, 'W': -4}, 
        'WWPSC': {'S': -2, 'W': -2}, 'WWPSA': {'S': -4, 'W': -4},
        'WSPWC': {'S': -2, 'W': -2}, 'WSPWA': {'S': -5, 'W': -5},
        'SWPWC': {'S': -2, 'W': -2}, 'SWPWA': {'S': -4, 'W': -4},
        'WSPSC': {'S': -2, 'W': -2}, 'WSPSA': {'S': -4, 'W': -4},
        'SWPSC': {'S': -2, 'W': -2}, 'SWPSA': {'S': -4, 'W': -4},
        'SSPWC': {'S': -4, 'W': -4}, 'SSPWA': {'S': -5, 'W': -5},
        'SSPSC': {'S': -2, 'W': -2}, 'SSPSA': {'S': -4, 'W': -4},

        # 3.2 B
        'WWNWC': {'S': -2, 'W': -2}, 'WWNWA': {'S': -4, 'W': -4}, 
        'WWNSC': {'S': -2, 'W': -4}, 'WWNSA': {'S': -5, 'W': -4},
        'WSNWC': {'S': -4, 'W': -2}, 'WSNWA': {'S': -5, 'W': -5},
        'SWNWC': {'S': -5, 'W': -2}, 'SWNWA': {'S': -4, 'W': -4},
        'WSNSC': {'S': -2, 'W': -2}, 'WSNSA': {'S': -4, 'W': -5},
        'SWNSC': {'S': -2, 'W': -2}, 'SWNSA': {'S': -5, 'W': -4},
        'SSNWC': {'S': -4, 'W': -4}, 'SSNWA': {'S': -5, 'W': -5},
        'SSNSC': {'S': -2, 'W': -2}, 'SSNSA': {'S': -4, 'W': -4},

        # 3.3 A
        'WPWWC': {'S': -2, 'W': -2}, 'WPWWA': {'S': -3, 'W': -3}, 
        'WPWSC': {'S': -2, 'W': -2}, 'WPWSA': {'S': -3, 'W': -3},
        'WPSWC': {'S': -2, 'W': -2}, 'WPSWA': {'S': -5, 'W': -5},
        'SPWWC': {'S': -2, 'W': -2}, 'SPWWA': {'S': -3, 'W': -3},
        'WPSSC': {'S': -2, 'W': -2}, 'WPSSA': {'S': -3, 'W': -3},
        'SPWSC': {'S': -2, 'W': -2}, 'SPWSA': {'S': -3, 'W': -3},
        'SPSWC': {'S': -3, 'W': -3}, 'SPSWA': {'S': -5, 'W': -5},
        'SPSSC': {'S': -2, 'W': -2}, 'SPSSA': {'S': -3, 'W': -3},

        # 3.3 B
        'WNWWC': {'S': -2, 'W': -2}, 'WNWWA': {'S': -3, 'W': -3},
        'WNWSC': {'S': -2, 'W': -3}, 'WNWSA': {'S': -5, 'W': -3},
        'WNSWC': {'S': -3, 'W': -2}, 'WNSWA': {'S': -5, 'W': -5},
        'SNWWC': {'S': -5, 'W': -2}, 'SNWWA': {'S': -3, 'W': -3},
        'WNSSC': {'S': -2, 'W': -2}, 'WNSSA': {'S': -3, 'W': -5},
        'SNWSC': {'S': -2, 'W': -2}, 'SNWSA': {'S': -5, 'W': -3},
        'SNSWC': {'S': -3, 'W': -3}, 'SNSWA': {'S': -5, 'W': -5},
        'SNSSC': {'S': -2, 'W': -2}, 'SNSSA': {'S': -3, 'W': -3},

        # 3.4 A
        'PWWWC': {'S': -2, 'W': -2}, 'PWWWA': {'S': -3, 'W': -3},
        'PWWSC': {'S': -2, 'W': -2}, 'PWWSA': {'S': -3, 'W': -3},
        'PWSWC': {'S': -2, 'W': -2}, 'PWSWA': {'S': -4, 'W': -4},
        'PSWWC': {'S': -2, 'W': -2}, 'PSWWA': {'S': -3, 'W': -3},
        'PWSSC': {'S': -2, 'W': -2}, 'PWSSA': {'S': -3, 'W': -3},
        'PSWSC': {'S': -2, 'W': -2}, 'PSWSA': {'S': -3, 'W': -3},
        'PSSWC': {'S': -3, 'W': -3}, 'PSSWA': {'S': -4, 'W': -4},
        'PSSSC': {'S': -2, 'W': -2}, 'PSSSA': {'S': -3, 'W': -3},

        # 3.4 B
        'NWWWC': {'S': -2, 'W': -2}, 'NWWWA': {'S': -3, 'W': -3},
        'NWWSC': {'S': -2, 'W': -3}, 'NWWSA': {'S': -4, 'W': -3},
        'NWSWC': {'S': -3, 'W': -2}, 'NWSWA': {'S': -4, 'W': -4},
        'NSWWC': {'S': -4, 'W': -2}, 'NSWWA': {'S': -3, 'W': -3},
        'NWSSC': {'S': -2, 'W': -2}, 'NWSSA': {'S': -3, 'W': -4},
        'NSWSC': {'S': -2, 'W': -2}, 'NSWSA': {'S': -4, 'W': -3},
        'NSSWC': {'S': -3, 'W': -3}, 'NSSWA': {'S': -4, 'W': -4},
        'NSSSC': {'S': -2, 'W': -2}, 'NSSSA': {'S': -3, 'W': -3},

    },

    frozenset(('G', 'T')) : {

        # 4.1 A
        'WWWPG': {'S': -4, 'W': -4}, 'WWWPT': {'S': -3, 'W': -3}, 
        'WWSPG': {'S': -3, 'W': -3}, 'WWSPT': {'S': -4, 'W': -4},
        'WSWPG': {'S': -4, 'W': -4}, 'WSWPT': {'S': -3, 'W': -3},
        'SWWPG': {'S': -4, 'W': -4}, 'SWWPT': {'S': -3, 'W': -3},
        'WSSPG': {'S': -4, 'W': -4}, 'WSSPT': {'S': -3, 'W': -3},
        'SWSPG': {'S': -3, 'W': -3}, 'SWSPT': {'S': -4, 'W': -4},
        'SSWPG': {'S': -4, 'W': -4}, 'SSWPT': {'S': -4, 'W': -4},
        'SSSPG': {'S': -4, 'W': -4}, 'SSSPT': {'S': -3, 'W': -3},

        # 4.1 B
        'WWWNG': {'S': -4, 'W': -4}, 'WWWNT': {'S': -3, 'W': -3}, 
        'WWSNG': {'S': -3, 'W': -5}, 'WWSNT': {'S': -4, 'W': -4},
        'WSWNG': {'S': -4, 'W': -5}, 'WSWNT': {'S': -3, 'W': -3},
        'SWWNG': {'S': -5, 'W': -4}, 'SWWNT': {'S': -3, 'W': -3},
        'WSSNG': {'S': -3, 'W': -4}, 'WSSNT': {'S': -3, 'W': -5},
        'SWSNG': {'S': -3, 'W': -5}, 'SWSNT': {'S': -3, 'W': -4},
        'SSWNG': {'S': -4, 'W': -5}, 'SSWNT': {'S': -4, 'W': -3},
        'SSSNG': {'S': -4, 'W': -4}, 'SSSNT': {'S': -3, 'W': -3},

        # 4.2 A
        'WWPWG': {'S': -2, 'W': -2}, 'WWPWT': {'S': -4, 'W': -4}, 
        'WWPSG': {'S': -2, 'W': -2}, 'WWPST': {'S': -4, 'W': -4},
        'WSPWG': {'S': -2, 'W': -2}, 'WSPWT': {'S': -5, 'W': -5},
        'SWPWG': {'S': -2, 'W': -2}, 'SWPWT': {'S': -4, 'W': -4},
        'WSPSG': {'S': -2, 'W': -2}, 'WSPST': {'S': -4, 'W': -4},
        'SWPSG': {'S': -2, 'W': -2}, 'SWPST': {'S': -4, 'W': -4},
        'SSPWG': {'S': -5, 'W': -5}, 'SSPWT': {'S': -4, 'W': -4},
        'SSPSG': {'S': -2, 'W': -2}, 'SSPST': {'S': -4, 'W': -4},

        # 4.2 B
        'WWNWG': {'S': -2, 'W': -2}, 'WWNWT': {'S': -4, 'W': -4}, 
        'WWNSG': {'S': -2, 'W': -5}, 'WWNST': {'S': -4, 'W': -4},
        'WSNWG': {'S': -4, 'W': -2}, 'WSNWT': {'S': -5, 'W': -5},
        'SWNWG': {'S': -5, 'W': -2}, 'SWNWT': {'S': -4, 'W': -4},
        'WSNSG': {'S': -2, 'W': -2}, 'WSNST': {'S': -4, 'W': -5},
        'SWNSG': {'S': -2, 'W': -2}, 'SWNST': {'S': -5, 'W': -4},
        'SSNWG': {'S': -5, 'W': -5}, 'SSNWT': {'S': -4, 'W': -4},
        'SSNSG': {'S': -2, 'W': -2}, 'SSNST': {'S': -4, 'W': -4},

        # 4.3 A
        'WPWWG': {'S': -2, 'W': -2}, 'WPWWT': {'S': -3, 'W': -3}, 
        'WPWSG': {'S': -2, 'W': -2}, 'WPWST': {'S': -3, 'W': -3},
        'WPSWG': {'S': -2, 'W': -2}, 'WPSWT': {'S': -5, 'W': -5},
        'SPWWG': {'S': -2, 'W': -2}, 'SPWWT': {'S': -3, 'W': -3},
        'WPSSG': {'S': -2, 'W': -2}, 'WPSST': {'S': -3, 'W': -3},
        'SPWSG': {'S': -2, 'W': -2}, 'SPWST': {'S': -3, 'W': -3},
        'SPSWG': {'S': -5, 'W': -5}, 'SPSWT': {'S': -3, 'W': -3},
        'SPSSG': {'S': -2, 'W': -2}, 'SPSST': {'S': -3, 'W': -3},

        # 4.3 B
        'WNWWG': {'S': -2, 'W': -2}, 'WNWWT': {'S': -3, 'W': -3},
        'WNWSG': {'S': -2, 'W': -5}, 'WNWST': {'S': -3, 'W': -3},
        'WNSWG': {'S': -3, 'W': -2}, 'WNSWT': {'S': -5, 'W': -5},
        'SNWWG': {'S': -5, 'W': -2}, 'SNWWT': {'S': -3, 'W': -3},
        'WNSSG': {'S': -2, 'W': -2}, 'WNSST': {'S': -3, 'W': -5},
        'SNWSG': {'S': -2, 'W': -2}, 'SNWST': {'S': -5, 'W': -3},
        'SNSWG': {'S': -5, 'W': -5}, 'SNSWT': {'S': -3, 'W': -3},
        'SNSSG': {'S': -2, 'W': -2}, 'SNSST': {'S': -3, 'W': -3},

        # 4.4 A
        'PWWWG': {'S': -2, 'W': -2}, 'PWWWT': {'S': -3, 'W': -3},
        'PWWSG': {'S': -2, 'W': -2}, 'PWWST': {'S': -3, 'W': -3},
        'PWSWG': {'S': -2, 'W': -2}, 'PWSWT': {'S': -4, 'W': -4},
        'PSWWG': {'S': -2, 'W': -2}, 'PSWWT': {'S': -3, 'W': -3},
        'PWSSG': {'S': -2, 'W': -2}, 'PWSST': {'S': -3, 'W': -3},
        'PSWSG': {'S': -2, 'W': -2}, 'PSWST': {'S': -3, 'W': -3},
        'PSSWG': {'S': -4, 'W': -4}, 'PSSWT': {'S': -3, 'W': -3},
        'PSSSG': {'S': -2, 'W': -2}, 'PSSST': {'S': -3, 'W': -3},

        # 4.4 B
        'NWWWG': {'S': -2, 'W': -2}, 'NWWWT': {'S': -3, 'W': -3},
        'NWWSG': {'S': -2, 'W': -4}, 'NWWST': {'S': -3, 'W': -3},
        'NWSWG': {'S': -3, 'W': -2}, 'NWSWT': {'S': -4, 'W': -4},
        'NSWWG': {'S': -4, 'W': -2}, 'NSWWT': {'S': -3, 'W': -3},
        'NWSSG': {'S': -2, 'W': -2}, 'NWSST': {'S': -3, 'W': -4},
        'NSWSG': {'S': -2, 'W': -2}, 'NSWST': {'S': -4, 'W': -3},
        'NSSWG': {'S': -4, 'W': -4}, 'NSSWT': {'S': -3, 'W': -3},
        'NSSSG': {'S': -2, 'W': -2}, 'NSSST': {'S': -3, 'W': -3},

    },

    frozenset(('G', 'A')) : {

        # 5.1 A
        'WWWPG': {'S': -3, 'W': -3}, 'WWWPA': {'S': -4, 'W': -4}, 
        'WWSPG': {'S': -3, 'W': -3}, 'WWSPA': {'S': -4, 'W': -4},
        'WSWPG': {'S': -4, 'W': -4}, 'WSWPA': {'S': -3, 'W': -3},
        'SWWPG': {'S': -3, 'W': -3}, 'SWWPA': {'S': -4, 'W': -4},
        'WSSPG': {'S': -3, 'W': -3}, 'WSSPA': {'S': -4, 'W': -4},
        'SWSPG': {'S': -3, 'W': -3}, 'SWSPA': {'S': -4, 'W': -4},
        'SSWPG': {'S': -4, 'W': -4}, 'SSWPA': {'S': -3, 'W': -3},
        'SSSPG': {'S': -3, 'W': -3}, 'SSSPA': {'S': -4, 'W': -4},

        # 5.1 B
        'WWWNG': {'S': -3, 'W': -3}, 'WWWNA': {'S': -4, 'W': -4}, 
        'WWSNG': {'S': -3, 'W': -4}, 'WWSNA': {'S': -5, 'W': -4},
        'WSWNG': {'S': -4, 'W': -3}, 'WSWNA': {'S': -5, 'W': -3},
        'SWWNG': {'S': -5, 'W': -3}, 'SWWNA': {'S': -4, 'W': -3},
        'WSSNG': {'S': -3, 'W': -3}, 'WSSNA': {'S': -4, 'W': -5},
        'SWSNG': {'S': -3, 'W': -3}, 'SWSNA': {'S': -5, 'W': -4},
        'SSWNG': {'S': -4, 'W': -4}, 'SSWNA': {'S': -5, 'W': -3},
        'SSSNG': {'S': -3, 'W': -3}, 'SSSNA': {'S': -4, 'W': -4},

        # 5.2 A
        'WWPWG': {'S': -2, 'W': -2}, 'WWPWA': {'S': -4, 'W': -4}, 
        'WWPSG': {'S': -2, 'W': -2}, 'WWPSA': {'S': -4, 'W': -4},
        'WSPWG': {'S': -2, 'W': -2}, 'WSPWA': {'S': -5, 'W': -5},
        'SWPWG': {'S': -2, 'W': -2}, 'SWPWA': {'S': -4, 'W': -4},
        'WSPSG': {'S': -2, 'W': -2}, 'WSPSA': {'S': -4, 'W': -4},
        'SWPSG': {'S': -2, 'W': -2}, 'SWPSA': {'S': -4, 'W': -4},
        'SSPWG': {'S': -4, 'W': -4}, 'SSPWA': {'S': -5, 'W': -5},
        'SSPSG': {'S': -2, 'W': -2}, 'SSPSA': {'S': -4, 'W': -4},

        # 5.2 B
        'WWNWG': {'S': -2, 'W': -2}, 'WWNWA': {'S': -4, 'W': -4}, 
        'WWNSG': {'S': -2, 'W': -4}, 'WWNSA': {'S': -5, 'W': -4},
        'WSNWG': {'S': -4, 'W': -2}, 'WSNWA': {'S': -5, 'W': -5},
        'SWNWG': {'S': -5, 'W': -2}, 'SWNWA': {'S': -4, 'W': -4},
        'WSNSG': {'S': -2, 'W': -2}, 'WSNSA': {'S': -4, 'W': -5},
        'SWNSG': {'S': -2, 'W': -2}, 'SWNSA': {'S': -5, 'W': -4},
        'SSNWG': {'S': -4, 'W': -4}, 'SSNWA': {'S': -5, 'W': -5},
        'SSNSG': {'S': -2, 'W': -2}, 'SSNSA': {'S': -4, 'W': -4},

        # 5.3 A
        'WPWWG': {'S': -2, 'W': -2}, 'WPWWA': {'S': -3, 'W': -3}, 
        'WPWSG': {'S': -2, 'W': -2}, 'WPWSA': {'S': -3, 'W': -3},
        'WPSWG': {'S': -2, 'W': -2}, 'WPSWA': {'S': -5, 'W': -5},
        'SPWWG': {'S': -2, 'W': -2}, 'SPWWA': {'S': -3, 'W': -3},
        'WPSSG': {'S': -2, 'W': -2}, 'WPSSA': {'S': -3, 'W': -3},
        'SPWSG': {'S': -2, 'W': -2}, 'SPWSA': {'S': -3, 'W': -3},
        'SPSWG': {'S': -3, 'W': -3}, 'SPSWA': {'S': -5, 'W': -5},
        'SPSSG': {'S': -2, 'W': -2}, 'SPSSA': {'S': -3, 'W': -3},

        # 5.3 B
        'WNWWG': {'S': -2, 'W': -2}, 'WNWWA': {'S': -3, 'W': -3},
        'WNWSG': {'S': -2, 'W': -3}, 'WNWSA': {'S': -5, 'W': -3},
        'WNSWG': {'S': -3, 'W': -2}, 'WNSWA': {'S': -5, 'W': -5},
        'SNWWG': {'S': -5, 'W': -2}, 'SNWWA': {'S': -3, 'W': -3},
        'WNSSG': {'S': -2, 'W': -2}, 'WNSSA': {'S': -3, 'W': -5},
        'SNWSG': {'S': -2, 'W': -2}, 'SNWSA': {'S': -5, 'W': -3},
        'SNSWG': {'S': -3, 'W': -3}, 'SNSWA': {'S': -5, 'W': -5},
        'SNSSG': {'S': -2, 'W': -2}, 'SNSSA': {'S': -3, 'W': -3},

        # 5.4 A
        'PWWWG': {'S': -2, 'W': -2}, 'PWWWA': {'S': -3, 'W': -3},
        'PWWSG': {'S': -2, 'W': -2}, 'PWWSA': {'S': -3, 'W': -3},
        'PWSWG': {'S': -2, 'W': -2}, 'PWSWA': {'S': -4, 'W': -4},
        'PSWWG': {'S': -2, 'W': -2}, 'PSWWA': {'S': -3, 'W': -3},
        'PWSSG': {'S': -2, 'W': -2}, 'PWSSA': {'S': -3, 'W': -3},
        'PSWSG': {'S': -2, 'W': -2}, 'PSWSA': {'S': -3, 'W': -3},
        'PSSWG': {'S': -3, 'W': -3}, 'PSSWA': {'S': -4, 'W': -4},
        'PSSSG': {'S': -2, 'W': -2}, 'PSSSA': {'S': -3, 'W': -3},

        # 5.4 B
        'NWWWG': {'S': -2, 'W': -2}, 'NWWWA': {'S': -3, 'W': -3},
        'NWWSG': {'S': -2, 'W': -3}, 'NWWSA': {'S': -4, 'W': -3},
        'NWSWG': {'S': -3, 'W': -2}, 'NWSWA': {'S': -4, 'W': -4},
        'NSWWG': {'S': -4, 'W': -2}, 'NSWWA': {'S': -3, 'W': -3},
        'NWSSG': {'S': -2, 'W': -2}, 'NWSSA': {'S': -3, 'W': -4},
        'NSWSG': {'S': -2, 'W': -2}, 'NSWSA': {'S': -4, 'W': -3},
        'NSSWG': {'S': -3, 'W': -3}, 'NSSWA': {'S': -4, 'W': -4},
        'NSSSG': {'S': -2, 'W': -2}, 'NSSSA': {'S': -3, 'W': -3},

    },

    frozenset(('A', 'T')) : {

        # 6.1 A
        'WWWPT': {'S': -3, 'W': -3}, 'WWWPA': {'S': -4, 'W': -4}, 
        'WWSPT': {'S': -4, 'W': -4}, 'WWSPA': {'S': -5, 'W': -5},
        'WSWPT': {'S': -3, 'W': -3}, 'WSWPA': {'S': -5, 'W': -5},
        'SWWPT': {'S': -3, 'W': -3}, 'SWWPA': {'S': -4, 'W': -4},
        'WSSPT': {'S': -3, 'W': -3}, 'WSSPA': {'S': -4, 'W': -4},
        'SWSPT': {'S': -3, 'W': -3}, 'SWSPA': {'S': -5, 'W': -5},
        'SSWPT': {'S': -4, 'W': -4}, 'SSWPA': {'S': -5, 'W': -5},
        'SSSPT': {'S': -3, 'W': -3}, 'SSSPA': {'S': -4, 'W': -4},

        # 6.1 B
        'WWWNT': {'S': -3, 'W': -3}, 'WWWNA': {'S': -4, 'W': -4}, 
        'WWSNT': {'S': -3, 'W': -4}, 'WWSNA': {'S': -3, 'W': -4},
        'WSWNT': {'S': -4, 'W': -3}, 'WSWNA': {'S': -4, 'W': -3},
        'SWWNT': {'S': -3, 'W': -3}, 'SWWNA': {'S': -4, 'W': -4},
        'WSSNT': {'S': -3, 'W': -3}, 'WSSNA': {'S': -4, 'W': -4},
        'SWSNT': {'S': -3, 'W': -4}, 'SWSNA': {'S': -3, 'W': -4},
        'SSWNT': {'S': -4, 'W': -3}, 'SSWNA': {'S': -4, 'W': -3},
        'SSSNT': {'S': -3, 'W': -3}, 'SSSNA': {'S': -4, 'W': -4},

        # 6.2 A
        'WWPWT': {'S': -4, 'W': -4}, 'WWPWA': {'S': -5, 'W': -5}, 
        'WWPST': {'S': -4, 'W': -4}, 'WWPSA': {'S': -5, 'W': -5},
        'WSPWT': {'S': -2, 'W': -2}, 'WSPWA': {'S': -5, 'W': -5},
        'SWPWT': {'S': -2, 'W': -2}, 'SWPWA': {'S': -4, 'W': -4},
        'WSPST': {'S': -2, 'W': -2}, 'WSPSA': {'S': -4, 'W': -4},
        'SWPST': {'S': -2, 'W': -2}, 'SWPSA': {'S': -5, 'W': -5},
        'SSPWT': {'S': -4, 'W': -4}, 'SSPWA': {'S': -5, 'W': -5},
        'SSPST': {'S': -4, 'W': -4}, 'SSPSA': {'S': -5, 'W': -5},

        # 6.2 B
        'WWNWT': {'S': -4, 'W': -4}, 'WWNWA': {'S': -5, 'W': -5}, 
        'WWNST': {'S': -4, 'W': -4}, 'WWNSA': {'S': -5, 'W': -5},
        'WSNWT': {'S': -4, 'W': -5}, 'WSNWA': {'S': -4, 'W': -5},
        'SWNWT': {'S': -5, 'W': -4}, 'SWNWA': {'S': -5, 'W': -4},
        'WSNST': {'S': -4, 'W': -5}, 'WSNSA': {'S': -4, 'W': -5},
        'SWNST': {'S': -5, 'W': -4}, 'SWNSA': {'S': -5, 'W': -4},
        'SSNWT': {'S': -4, 'W': -4}, 'SSNWA': {'S': -5, 'W': -5},
        'SSNST': {'S': -4, 'W': -4}, 'SSNSA': {'S': -5, 'W': -5},

        # 6.3 A
        'WPWWT': {'S': -3, 'W': -3}, 'WPWWA': {'S': -5, 'W': -5}, 
        'WPWST': {'S': -3, 'W': -3}, 'WPWSA': {'S': -5, 'W': -5},
        'WPSWT': {'S': -2, 'W': -2}, 'WPSWA': {'S': -5, 'W': -5},
        'SPWWT': {'S': -2, 'W': -2}, 'SPWWA': {'S': -3, 'W': -3},
        'WPSST': {'S': -2, 'W': -2}, 'WPSSA': {'S': -3, 'W': -3},
        'SPWST': {'S': -2, 'W': -2}, 'SPWSA': {'S': -5, 'W': -5},
        'SPSWT': {'S': -3, 'W': -3}, 'SPSWA': {'S': -5, 'W': -5},
        'SPSST': {'S': -3, 'W': -3}, 'SPSSA': {'S': -5, 'W': -5},

        # 6.3 B
        'WNWWT': {'S': -3, 'W': -3}, 'WNWWA': {'S': -5, 'W': -5},
        'WNWST': {'S': -3, 'W': -3}, 'WNWSA': {'S': -5, 'W': -5},
        'WNSWT': {'S': -3, 'W': -5}, 'WNSWA': {'S': -3, 'W': -5},
        'SNWWT': {'S': -5, 'W': -3}, 'SNWWA': {'S': -5, 'W': -3},
        'WNSST': {'S': -3, 'W': -5}, 'WNSSA': {'S': -3, 'W': -5},
        'SNWST': {'S': -5, 'W': -3}, 'SNWSA': {'S': -5, 'W': -3},
        'SNSWT': {'S': -3, 'W': -3}, 'SNSWA': {'S': -5, 'W': -5},
        'SNSST': {'S': -3, 'W': -3}, 'SNSSA': {'S': -5, 'W': -5},

        # 6.4 A
        'PWWWT': {'S': -3, 'W': -3}, 'PWWWA': {'S': -4, 'W': -4},
        'PWWST': {'S': -3, 'W': -3}, 'PWWSA': {'S': -4, 'W': -4},
        'PWSWT': {'S': -2, 'W': -2}, 'PWSWA': {'S': -4, 'W': -4},
        'PSWWT': {'S': -2, 'W': -2}, 'PSWWA': {'S': -3, 'W': -3},
        'PWSST': {'S': -2, 'W': -2}, 'PWSSA': {'S': -3, 'W': -3},
        'PSWST': {'S': -2, 'W': -2}, 'PSWSA': {'S': -4, 'W': -4},
        'PSSWT': {'S': -3, 'W': -3}, 'PSSWA': {'S': -4, 'W': -4},
        'PSSST': {'S': -3, 'W': -3}, 'PSSSA': {'S': -4, 'W': -4},

        # 6.4 B
        'NWWWT': {'S': -3, 'W': -3}, 'NWWWA': {'S': -4, 'W': -4},
        'NWWST': {'S': -3, 'W': -3}, 'NWWSA': {'S': -4, 'W': -4},
        'NWSWT': {'S': -3, 'W': -4}, 'NWSWA': {'S': -3, 'W': -4},
        'NSWWT': {'S': -4, 'W': -3}, 'NSWWA': {'S': -4, 'W': -3},
        'NWSST': {'S': -3, 'W': -4}, 'NWSSA': {'S': -3, 'W': -4},
        'NSWST': {'S': -4, 'W': -3}, 'NSWSA': {'S': -4, 'W': -3},
        'NSSWT': {'S': -3, 'W': -3}, 'NSSWA': {'S': -4, 'W': -4},
        'NSSST': {'S': -3, 'W': -3}, 'NSSSA': {'S': -4, 'W': -4},

        }}

def differences_at_end(pair, snp_position):
    """
    If snp_position is 'last', return the number of nucleotide
    differences in the last 4 nucleotides. If snp_position is
    'first', return the number of differences in the first 4
    nucleotides.

    Args:
        pair: A 2-tuple of AMAS primers
        snp_position: Relative SNP position. Either 'first' or
            'last'.
    
    Returns:
        The hamming distance between the first or last 4 nucleotides.
    
    Raise:
        ValueError if either primer is shorter than 4 nucleotides.
    """

    seq1 = str(pair[0].sequence)
    seq2 = str(pair[1].sequence)

    if len(seq1) < 4 or len(seq2) < 4:
        raise ValueError('Primers must be greater than 4 nucleotides.')

    if snp_position == 'first':
        return hamming(seq1[:4], seq2[:4])
    else:
        return hamming(seq1[-4:], seq2[-4:])

def preserve_best_and_substitute(pairs, snp_position):
    """
    Corresponds to One SNP module, Two SNP module, and Three SNP module
    from STARP F primer design_clarified.docx

    (1) Find the number of SNPs at the end of the shortest primer pair.
        This value influences the cutoff used when checking if the total
        SNP number between each pair is greater than 'cutoff'.

    (2) The temperature range and desired average temperature is chosen
        based on if there are any 'high snp pairs' from part 1. If there
        are none, then the substitute flag is set to true to signal that
        the selected primer pair will undergo substitutions.

    If no primers fit these conditions or the list of pairs is empty,
    None is returned.

    Args:
        pairs: A list of 2-tuples of AMAS pairs.
        snp_position: Position of the SNP relative to the AMAS primers.
            Should be either 'first' or 'last'.

    Returns:
        A 2-tuple representing the best AMAS pair.
        None if no primers have the described qualities.
    """

    # If the list is empty there is no good primer.
    if not pairs:
        return None

    # Flag for checking if substitutions need to be made on the selected
    # AMAS primer.
    substitute = False

    # Sort pairs by their length, shortest to longest.
    pairs = sorted(pairs, key=len)

    # Calculate SNP numbers at the ends of the first pair.
    snp_number = differences_at_end(pairs[0], snp_position)

    # All pairs should have 1 difference at the SNP location.
    if snp_number == 1:
        total_snp_number_cutoff = 4
    elif snp_number == 2:
        total_snp_number_cutoff = 3
    else:
        # All primers will be selected in this case.
        total_snp_number_cutoff = 0

    # Get pairs with >= 'total_snp_number_cutoff' SNPs
    high_snp_pairs = list(filter(
        lambda pair: hamming(pair[0].sequence, pair[1].sequence) >= total_snp_number_cutoff,
        pairs
    ))

    if high_snp_pairs:
        if snp_number > 2:
            low = 52
            high = 58
            average = 58
        else:
            # If there are some, try find the pair with Tm between 53
            # and 60 C and average is closest to 53 C.
            low = 53
            high = 60
            average = 53
    else:
        substitute = True
        low = 54
        high = 58
        average = 58

    # Keep pairs with Tm between 'low' and 'high'
    if high_snp_pairs:
        acceptable_pairs = list(filter(
            lambda pair: low <= pair[0].tm <= high and low <= pair[1].tm <= high,
            high_snp_pairs
        ))
    else:
        acceptable_pairs = list(filter(
            lambda pair: low <= pair[0].tm <= high and low <= pair[1].tm <= high,
            pairs
        ))

    # Sort pairs by average Tm closest to 'average' C.
    sorted_pairs = sorted(
        acceptable_pairs,
        key=lambda pair: abs(((pair[0].tm + pair[1].tm) / 2) - average)
    )

    # If there are acceptable primer pairs, select the first one from
    # the sorted list since it is closest to the desired temperature.
    # If there are no acceptable primer pairs, then select the primers
    # between low and high, closest to average even if they are not
    # paired with each other. Return None if still no acceptable primer
    # pair can be found.
    if sorted_pairs:
        best_pair = sorted_pairs[0]
    else:
        allele1_primers = sorted(
            [pair[0] for pair in pairs],
            key=lambda primer: abs(primer.tm - average)
        )
        allele2_primers = sorted(
            [pair[1] for pair in pairs],
            key=lambda primer: abs(primer.tm - average)
        )

        if low <= allele1_primers[0].tm <= high and low <= allele2_primers[0].tm <= high:
            best_pair = (allele1_primers[0], allele2_primers[0])
        else:
            return None

    # If the flag is on, the bases of the best pair need to be
    # substituted.
    if substitute:
        best_pair[0].sequence, best_pair[1].sequence = substitute_bases((best_pair[0].sequence, best_pair[1].sequence), snp_position=snp_position)
    
    return best_pair

def amas_pair_filter(amas_pairs, snp_position):
    """
    Removes primers pairs when one of the constituent primers has:
    1. >= 10 contiguous G/C or >= 12 contiguous A/T
    2. Mononucleotide repeat of length 8+
    3. If snp_position == 'first':
           has mononucleotide repeat of 4 A/Ts or 5 G/Cs in first 6 bases
       Else if snp_position == 'last':
           has mononucleotide repeat of 4 A/Ts or 5 G/Cs in last 6 bases
    4. Dinucleotide repeat of length 6+
    5. GC > 0.80 or GC < 0.20

    Args:
        amas_pairs: A list of 2-tuples of AmasPrimers.
    """
    def filter_func(pair):
        for primer in pair:
            seq = Sequence(primer.sequence)
            if (seq.has_contig_gc_at(10, 12)
                    or seq.has_mononucleotide_repeat(8, 8)
                    or seq.has_dinucleotide_repeat(6)
                    or seq.gc < 0.20
                    or seq.gc > 0.80):
                return False
            
            if snp_position == 'last' and seq[-6:].has_mononucleotide_repeat(5, 4):
                return False

            if snp_position == 'first' and seq[:6].has_mononucleotide_repeat(5, 4):
                return False

        return True

    return list(filter(filter_func, amas_pairs))

def substitute(seq, idx):
    """
    Substitutes the given index of the allele according to the
    mapping nucleotide_sub. This is necessary because strings
    and Sequence objects are immutable.

    Args:
        seq: the Sequence object to modify
        idx: an integer or tuple of the indices to change.

    Returns:
        The modified allele as a Sequence object.
    """
    if isinstance(idx, int):
        idx = (idx,)  # Change int to 1-tuple.

    seq = str(seq)

    nucleotide_sub = {'A' : 'C', 'T' : 'C', 'G' : 'A', 'C' : 'T'}
    seq = list(seq)

    for i in idx:
        seq[i] = nucleotide_sub[seq[i]]

    return Sequence(''.join(seq))

def seq_to_ambiguity_code(sequence: str):
    """ Converts 'C' and 'G' to 'S', and 'A'/'T' to 'W' as defined by
    http://www.reverse-complement.com/ambiguity.html


    Note that these are not being used to designate SNPs. In Dr. Long's
    instructions, many times he asks for a count of how many G/C bases
    are in a certain sequence, and this function reduces the sizes of
    the already massive dictionaries above.
    """
    return (sequence.replace('C', 'S').replace('G', 'S')
            .replace('A', 'W').replace('T', 'W'))

def generate_amas_for_substitution(allele1, allele2, position):
    """ Attempts to find the best upstream pair and the best
    downstream pair. These could be None.
    
    Args:
        allele1: The aligned first allele.
        allele2: The aligned second allele.
        position: The index around which to generate primers.
    """
    pairs = list(zip(generate_amas_upstream(allele1, 1, position, 16, 26),
                     generate_amas_upstream(allele2, 2, position, 16, 26)))

    # Remove pairs whose primers have undesirable characteristics.
    pairs = amas_pair_filter(pairs, snp_position='last')

    upstream_pair = preserve_best_and_substitute(pairs, snp_position='last')

    pairs = list(zip(generate_amas_downstream(allele1, 1, position, 16, 26),
                     generate_amas_downstream(allele2, 2, position, 16, 26)))

    pairs = amas_pair_filter(pairs, snp_position='first')
    
    downstream_pair = preserve_best_and_substitute(pairs, snp_position='first')

    return upstream_pair, downstream_pair

def generate_amas_for_indel(allele1, allele2, position):
    """
    
    An example works best.

    Say we have the alleles

        GTGG ACGCTCGAGGACTATAG--TCAGGAGAGGTGGGCATGG
        |||| |||||||||||||||||  |||||||||||||||||||
        GTGG ACGCTCGAGGACTATAGTCTCAGGAGAGGTGGGCATGG

    Then the upstream pairs that get generated are

        ACGCTCGAGGACTATAGT
        ||||||||||||||||||
        ACGCTCGAGGACTATAGT

        ACGCTCGAGGACTATAGTC
        |||||||||||||||||||
        ACGCTCGAGGACTATAGTC

        ACGCTCGAGGACTATAGTCA
        |||||||||||||||||||
        ACGCTCGAGGACTATAGTCT

        ACGCTCGAGGACTATAGTCAG
        |||||||||||||||||||
        ACGCTCGAGGACTATAGTCTC

        ACGCTCGAGGACTATAGTCAGG
        |||||||||||||||||||
        ACGCTCGAGGACTATAGTCTCA

        ...

        ACGCTCGAGGACTATAGTCAGGAGA
        |||||||||||||||||||    ||
        ACGCTCGAGGACTATAGTCTCAGGA

    Args:
        allele1: The aligned first allele.
        allele2: The aligned second allele.
        position: The index of the indel around which to generate
            primers.

    Returns:
        upstream_pair, downstream_pair

        upstream_pair: The best AMAS pair where the majority of the
            primer sequence is upstream of the SNP.
        downstream_pair: The best AMAS pair where the majority of the
            primer sequence is downstream from the SNP.
    """

    # Create upstream AMAS primers by moving the SNP back 15 positions
    # and creating downstream primers from that new position.
    pairs = zip(generate_amas_downstream(allele1, 1, position-15, 16, 26),
                     generate_amas_downstream(allele2, 2, position-15, 16, 26))

    # Remove pairs with the same nucleotide on the 3' end.
    pairs = list(filter(lambda pair: pair[0][-1] != pair[1][-1], pairs))

    # Filter pairs for undesirable characteristics.
    pairs = amas_pair_filter(pairs, snp_position='last')

    # Order the pairs based on
    # 1) The number of nucleotide differences in the last 4 bases.
    #    More differences is preferable.
    # 2) Length of the sequences in each pair. Shorter is preferable.
    upstream_pairs = sorted(pairs,
                            key=lambda pair: (
                                Sequence.hamming(pair[0][-4:], pair[1][-4:]),
                                1 / len(pair[0])),
                            reverse=True)

    # Go through the list, verifying the melting temperatures are good and
    # substituting bases if needed. Since these are already sorted, the first
    # primer that succeeds is the one sought after.
    upstream_pair = None
    for pair in upstream_pairs:
        upstream_pair = preserve_best_and_substitute([pair], snp_position='last')
        if upstream_pair:
            break

    # Create downstream AMAS primers by moving the SNP forward 15
    # positions and creating upstream primers from that position.
    pairs = zip(generate_amas_upstream(allele1, 1, position+15, 16, 26),
                generate_amas_upstream(allele2, 2, position+15, 16, 26))

    # Remove pairs with the same nucleotide on the 5' end.
    pairs = list(filter(lambda pair: pair[0][0] != pair[1][0], pairs))

    # Filter pairs for undesirable characteristics
    pairs = amas_pair_filter(pairs, snp_position='first')

    # Order the pairs based on
    # 1) The number of nucleotide differences in the first 4 bases.
    #    More differences is preferable.
    # 2) Length of the sequences in each pair. Shorter is preferable.
    downstream_pairs = sorted(pairs,
                              key=lambda pair: (
                                  Sequence.hamming(pair[0][:4], pair[1][:4]),
                                  1 / len(pair[0])),
                              reverse=True)

    # Go through the list, verifying the melting temperatures are good and
    # substituting bases if needed. Since these are already sorted, the first
    # primer that succeeds is the one sought after.
    downstream_pair = None
    for pair in downstream_pairs:
        downstream_pair = preserve_best_and_substitute([pair], snp_position='first')
        if downstream_pair:
            break

    return upstream_pair, downstream_pair

def generate_amas_upstream(allele, allele_num, pos, minimum, maximum):
    """ Returns AMAS primers upstream of the position using the 
    aligned sequence. 

    Args:
        allele: The aligned allele to create primers from. It is
            possible for this to contain dashes.
        allele_num: This is using either allele 1 or 2. This is
            necessary for the primer instantiation.
        pos: The position to start from.
        minimum: The minimum primer length.
        maximum: The maximum primer length, inclusive.
    """
    sliced = str(allele[:pos+1]).replace('-', '')
    return [AmasPrimer(sliced[0-size:],
                       allele_num,
                       span=(len(sliced)-size, len(sliced)))
            for size in range(minimum, maximum+1)
            if pos >= size]

def generate_amas_downstream(allele, allele_num, pos, minimum, maximum):
    """ Returns a list of AMAS primers downstream of the position using
    the aligned allele. Doing this downstream is slightly more
    complicated since we have to consider the positions of the dashes.
    """
    dashes_before_pos = str(allele[:pos]).count('-')
    dashes_after_pos = str(allele[pos+1:]).count('-')
    is_dash = 1 if allele[pos] == '-' else 0

    sliced = str(allele[pos:]).replace('-', '')
    return [AmasPrimer(sliced[:size],
                       allele_num,
                       span=(pos-dashes_before_pos,
                             pos-dashes_before_pos+size))
            for size in range(minimum, maximum+1)
            if size <= len(sliced)]

def substitute_bases(pair, snp_position='last'):
    """
    Substitutes bases according to Dr. Long's instructions.
    See 'docs/starp/STARP F primer design_clarified.docx'.

    Args:
        pair: A 2-tuple of Sequence objects or strings.
        snp_position: should be either 'first' or 'last', signifying
            if the snp is at the beginning or end of the pair.
    """
    if not pair:
        return None

    if not ((type(pair[0]) == Sequence and type(pair[1]) == Sequence)
            or (type(pair[0]) == str and type(pair[1]) == str)):
        raise ValueError('Pair must be a tuple of strings or Sequences')

    # All of the keys in the large dicts require 5+ characters.
    if len(pair[0]) < 5:
        return pair

    if snp_position == 'first':
        str_pair = (str(Sequence(pair[0]).rev_comp()), str(Sequence(pair[1]).rev_comp()))
    else:
        str_pair = (str(pair[0]), str(pair[1]))

    snp = Snp(f'.{len(str_pair[0])-1}{str_pair[0][-1]}>{str_pair[1][-1]}')

    # Record the SNPs in the 2nd, 3rd, 4th, or 5th position from the 3' end
    local_snps = TwoAlleles(f'>\n{str_pair[0][-5:-1]}\n>\n{str_pair[1][-5:-1]}').snps()

    if len(local_snps) == 0:
        new_amas1, new_amas2 = substitute_with_one_snp(str_pair, 'last')
        pair = (str(new_amas1), str(new_amas2))
    elif len(local_snps) == 1:
        new_amas1, new_amas2 = substitute_with_two_snps(str_pair, snp, 'last')
        pair = (str(new_amas1), str(new_amas2))

    pair = (Sequence(pair[0]), Sequence(pair[1]))

    # Reorient primers
    if snp_position == 'first':
        pair = (pair[0].rev_comp(), pair[1].rev_comp())

    return pair

def substitute_with_one_snp(pair, snp_position='last'):
    """
    Substitute the bases of the pair's sequences when the only SNP
    occurs at either the first or last nucleotide. This function probably
    should only be called by substitute_bases().

    This function assumes the primers in the pair are oriented on the
    plus strand.

    Notation for ambiguity codes comes from
    http://www.reverse-complement.com/ambiguity.html

    Ambiguity codes:
    G/C = S
    A/T = W

    Args:
        pair: A tuple of same-length strings or sequences.
        snp_position: Acceptable values are 'first' and 'last'.

    Returns:
        A 2-tuple of the pair with appropriate bases substituted
            of type Sequence.
    """

    # If the pair are downstream primers, as noted by snp_position
    # being 'first', we can reverse complement them and still use
    # the same instructions.
    if snp_position == 'first':
        pair = (str(pair[0].rev_comp()), str(pair[1].rev_comp()))
    else:
        pair = (str(pair[0]), str(pair[1]))

    if not pair[0][-1] != pair[1][-1]:
        raise ValueError('The sequences do not have a SNP in the last '
                         'position.')

    if pair[0][-4:-1] != pair[1][-4:-1]:
        raise ValueError('The sequences must be equal in the 2nd, 3rd, and '
                         '4th position from the 3\' end.')

    snp = Snp(f'.{len(pair[0])}{pair[0][-1]}>{pair[1][-1]}')

    # The only Snp between the two sequences is at the last index.
    code = seq_to_ambiguity_code(pair[0][-4:-1]) + pair[0][-1]
    idx_to_sub = sub_index_one_snp[snp.nucleotides][code]
    seq1 = substitute(pair[0], idx_to_sub)

    code = seq_to_ambiguity_code(pair[1][-4:-1]) + pair[1][-1]
    idx_to_sub = sub_index_one_snp[snp.nucleotides][code]
    seq2 = substitute(pair[1], idx_to_sub)

    # Orient the sequences back to their original orientation.
    if snp_position == 'first':
        seq1 = seq1.rev_comp()
        seq2 = seq2.rev_comp()

    return (seq1, seq2)

def substitute_with_two_snps(pair, snp, snp_position='last'):
    """
    Substitute the bases of the pair sequences when there are two SNPs
    between the sequences in the last five bases. This function probably
    should only be called by substitute_bases().

    This function assumes the primers in the pair are oriented on the
    plus strand.

    For reference, see
    docs/starp/STARP F primer design_clarified.docx

    Long's written instructions assume the SNP is placed at the end of
    the sequences. However, this function will also be called with a
    SNP at the beginning of the sequences. In this case, his
    verbal instructions were to perform the same operations but at the
    beginning of the sequences. To avoid making another very similar
    function, the sequences are reversed here then returned to their
    original orientation at the end of the function.
    """

    # If the pair are downstream primers, as noted by snp_position
    # being 'first', we can reverse complement them and still use
    # the same instructions.
    if snp_position == 'first':
        pair = (str(pair[0].rev_comp()), str(pair[1].rev_comp()))
    else:
        pair = (str(pair[0]), str(pair[1]))

    if not pair[0][-1] != pair[1][-1]:
        raise ValueError('The sequences do not have a SNP in the last '
                         'position.')

    # All SNPs between the two alleles.
    snps = TwoAlleles(f'>\n{pair[0]}\n>\n{pair[1]}').snps()

    # SNPs at 2nd, 3rd, 4th, or 5th position from the 3' end.
    snps = list(filter(lambda snp: 0 < len(pair[0])-snp.position-1 < 5, snps))

    if len(snps) != 1:
        raise ValueError('The sequences must have one SNP in the 2nd, 3rd, 4th, or '
                         '5th position from the 3\' end.')

    xsnp = snps[0]  # extra snp

    # Part of the codes created later concerns themselves about whether
    # the additional SNP is a [C/G] or [A/T] SNP, or the others. To simplify the possible
    # number of codes, [C/G] and [A/T] SNPs are designated as 'P' for paired
    # and all the rest are designated 'N'.
    if xsnp.nucleotides == {'C', 'G'} or xsnp.nucleotides == {'A', 'T'}:
        placeholder = 'P'
    else:
        placeholder = 'N'

    encoded_seq = pair[0][:xsnp.position] + placeholder + pair[0][xsnp.position+1:]
    code = seq_to_ambiguity_code(encoded_seq[-5:-1]) + pair[0][-1]
    idx_to_sub = sub_index_two_snps[snp.nucleotides][code][seq_to_ambiguity_code(xsnp.ref_nucleotide)]

    seq1 = substitute(pair[0], idx_to_sub)

    encoded_seq = pair[1][:xsnp.position] + placeholder + pair[1][xsnp.position+1:]
    code = seq_to_ambiguity_code(encoded_seq[-5:-1]) + pair[1][-1]
    idx_to_sub = sub_index_two_snps[snp.nucleotides][code][seq_to_ambiguity_code(xsnp.new_nucleotide)]

    seq2 = substitute(pair[1], idx_to_sub)

    # Orient the sequences back to their original orientation.
    if snp_position == 'first':
        seq1 = seq1.rev_comp()
        seq2 = seq2.rev_comp()

    return (seq1, seq2)
