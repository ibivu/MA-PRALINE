import regex as re

from pyparsing import OneOrMore, ZeroOrMore, Suppress, Or, Literal, Group, Dict
from pyparsing import Optional, Word

_SYMBOLS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
            'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
            'W', 'X', 'Y', 'Z']
_NUMBERS = '0123456789'
_SYMBOL_MATCH_ANY = 'x'

def _convert_integers(tokens):
    return int(tokens[0])

ParOpen = Suppress(Literal("("))
ParClose = Suppress(Literal(")"))
SqBrOpen = Suppress(Literal("["))
SqBrClose = Suppress(Literal("]"))
CuBrOpen = Suppress(Literal("{"))
CuBrClose = Suppress(Literal("}"))
Comma = Suppress(Literal(","))

AnyMatcher = Group(Literal(_SYMBOL_MATCH_ANY)).setResultsName("AnySymbol")
SymbolMatcher = Group(Or(Literal(symbol) for symbol in _SYMBOLS))
SymbolMatcher = SymbolMatcher.setResultsName("Symbol")
SymbolOrAnyMatcher = Or([SymbolMatcher, AnyMatcher])

NotInMatcher = Group(CuBrOpen + OneOrMore(SymbolMatcher) + CuBrClose)
NotInMatcher = NotInMatcher.setResultsName("NotInClass")
InMatcher = Group(SqBrOpen + OneOrMore(SymbolMatcher) + SqBrClose)
InMatcher = InMatcher.setResultsName("InClass")

PosInt = Word(_NUMBERS).setParseAction(_convert_integers)
TimesMatchers = Or([SymbolOrAnyMatcher, NotInMatcher, InMatcher])
NTimesQualifier = ParOpen + PosInt + ParClose
NTimesMatcher = Group(TimesMatchers + NTimesQualifier)
NTimesMatcher = NTimesMatcher.setResultsName("NTimes")
MinMaxTimesQualifier = ParOpen + PosInt + Comma + PosInt + ParClose
MinMaxTimesMatcher = Group(TimesMatchers + MinMaxTimesQualifier)
MinMaxTimesMatcher = MinMaxTimesMatcher.setResultsName("MinMaxTimes")

AllMatchers = [SymbolOrAnyMatcher, NotInMatcher, InMatcher, NTimesMatcher,
               MinMaxTimesMatcher]
PrositeElement = Or(AllMatchers)
Separator = Suppress(Literal("-"))

NTerminalMatcher = Optional(Group(Literal("<")).setResultsName("NTerminal"))
CTerminalMatcher = Optional(Group(Literal(">")).setResultsName("CTerminal"))

PrositePattern = NTerminalMatcher + PrositeElement + \
                 ZeroOrMore(Separator + PrositeElement) + CTerminalMatcher
PrositePattern = PrositePattern.setResultsName("Pattern")

def main():
    pattern = "<{ABC}(2)-x-x-x-F-[AB]-[CD](2,3)>"

    pat_tokens = PrositePattern.parseString(pattern)
    print pat_tokens
    for token in pat_tokens:
        print token.getName()

    print pattern_to_re(pattern)

def pattern_to_re(pattern):
    convertor = ParseTreeConvertor()
    pattern_tokens = PrositePattern.parseString(pattern)
    
    return re.compile(convertor(pattern_tokens))

class ParseTreeConvertor(object):
    def __init__(self):
            self.match_pos = 0

    def __call__(self, tree, set_capture_group=True):
        name = tree.getName()
        if name == "Pattern":
            for subtree in tree:
                return "".join([self(subtree) for subtree in tree])         
        elif name == "NTerminal":
            return "^"
        elif name == "CTerminal":
            return "$"
        elif name == "AnySymbol":
            return "(?P<spacer>[A-Z])"
        elif name == "Symbol":
            if set_capture_group:
                self.match_pos += 1
                fmt = "(?P<match_{0}>{1})"
                return fmt.format(self.match_pos, tree[0])
            else:
                return tree[0]
            return fmt.format(tree[0])
        elif name == "NotInClass":
            fmt = "(?P<spacer>[^{0}])"
            syms = [self(subtree, False) for subtree in tree]
            return fmt.format("".join(syms))
        elif name == "InClass":
            self.match_pos += 1
            fmt = "(?P<match_{0}>[{1}])"
            syms = [self(subtree, False) for subtree in tree]
            return fmt.format(self.match_pos, "".join(syms))
        elif name == "NTimes":
            fmt = "{0}{{{1}}}"
            return fmt.format(self(tree[0]), tree[1])
        elif name == "MinMaxTimes":
            fmt = "{0}{{{1},{2}}}"
            return fmt.format(self(tree[0]), tree[1], tree[2])

if __name__ == '__main__':
    main()