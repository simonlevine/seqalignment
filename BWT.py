# BWT.py

def rle(s):
    """Run Length Encoder
    Args: s, string to be encoded
    Returns: RLE(s)
    """

    def countRL(str, k):
        currcount = 0
        for char in str[k:]:
            if str[k] == char:
                currcount += 1
            else:
                break  # stop scanning to save effort

        return currcount

    for i in range(0, len(s)-1):

        if countRL(s, i) >= 2:

            rl = countRL(s, i)  # count of consecutive characters at s[i]
            #now to remove characters

            s = s[:i+2]+str(rl)+s[i+rl:]

    return s


def bwt_encode(s):
    """Burrows-Wheeler Transform
    Args: s, string, which must not contain '{' or '}'
    Returns: BWT(s), which contains '{' and '}'
    """
    s = "{" + s + "}"  # adding marker to start/end of text

    # generate an array of circularly-rotated strings, sorted,
    table = sorted(s[i:] + s[:i] for i in range(len(s)))

    #get last characters from each row.
    last_column = [row[-1:] for row in table]

    # get string of last characters in each column.
    return "".join(last_column)


def bwt_decode(bwt):
    """Inverse Burrows-Wheeler Transform
    Args: bwt, BWT'ed string, which should contain '{' and '}'
    Returns: reconstructed original string s, must not contains '{' or '}'
    """
    table = [""] * len(bwt)  # make an empty table
 
    for i in range(len(bwt)):
        table = sorted(bwt[i] + table[i]
                       for i in range(len(bwt)))  # Add a column of r
    # Find the correct row (ending in '}')

    s = [row for row in table if row.endswith("}")][0]

    return s.rstrip("}").strip("{")  # Get rid of start and end markers


def test_string(s):
    compressed = rle(s)
    bwt = bwt_encode(s)
    compressed_bwt = rle(bwt)
    reconstructed = bwt_decode(bwt)
    template = "{:25} ({:3d}) {}"
    print(template.format("original", len(s), s))
    print(template.format("bwt_enc(orig)", len(bwt), bwt))
    print(template.format("bwt_dec(bwt_enc(orig))", len(reconstructed), reconstructed))
    print(template.format("rle(orig)", len(compressed), compressed))
    print(template.format("rle(bwt_enc(orig))", len(compressed_bwt), compressed_bwt))
    print()
    print()


if __name__ == "__main__":
    # Add more of your own strings to explore for question (i)
    test_strings = ["WOOOOOHOOOOHOOOO!",
                    "scottytartanscottytartanscottytartanscottytartan"]
    for s in test_strings:
        test_string(s)
