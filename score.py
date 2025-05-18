def spscore(strs: list, d=3, e=1, m=1, mis=-2):
    nums = len(strs)
    score = 0
    for i in range(nums):
        for j in range(i + 1, nums):
            score += _spTwo(strs[i], strs[j])

    return score


def _spTwo(s1: str, s2: str, d=3, e=1, m=1, mis=-2):
    if len(s1) != len(s2):
        print(s1, s2)
        raise ValueError("seqs not aligned!")
    score_two = 0
    gap1 = 0
    for i in range(len(s1)):
        if s1[i] != "-" and s2[i] != "-":
            if gap1:
                gap1 = 0
            if s1[i] == s2[i]:
                score_two += m
            else:
                score_two += mis
        elif s1[i] != "-" or s2[i] != "-":
            if gap1 == 0:
                score_two -= d
                gap1 = 1
            else:
                score_two -= e
    return score_two