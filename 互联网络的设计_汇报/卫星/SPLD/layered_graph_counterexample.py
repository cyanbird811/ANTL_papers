# -*- coding: utf-8 -*-
"""
SPLD / MEGAREDUCE 可行性判定的松弛性验证

问题：论文 Algorithm 2 用「有向分层图上的单位容量最大流」判定
      (a) r 条边不相交路径  (b) 每条路径 <= L 跳
      但分层图里同一条物理边 (i,j) 在每一层都是一条独立的弧 (i,j,l)，
      单位容量只保证「分层弧不相交」，不保证「物理链路不相交」。

本脚本对比：
      mf   = 分层图 max-flow                     （论文判定用的量）
      true = 真实物理边不相交、且每条 <= L 跳 的路径条数（暴力枚举）
      若 mf > true，说明算法把不可行判成了可行 → 判定是原问题的松弛。
"""

import itertools
import random
from collections import deque


# ---------------------------------------------------------------- 原图
# 例 1：s-x-a-b-t 链 + 捷径 s-a, s-x
DEFAULT_EDGES = [('s', 'a'), ('s', 'x'), ('x', 'a'), ('a', 'b'), ('b', 't')]


# ------------------------------------------------- 论文 Algorithm 2 的做法
def build_layered(edges, L, src, dst):
    """§III-C 的分层图构造，返回弧集合（每条弧容量 1）"""
    V = set()
    for u, v in edges:
        V.update((u, v))
    layers = [set() for _ in range(L + 2)]      # 索引 1..L+1
    layers[1] = {src}
    for l in range(2, L + 1):
        layers[l] = V - {src}
    layers[L + 1] = {dst}

    arcs = []
    for l in range(1, L + 1):
        for u, v in edges:
            for a, b in ((u, v), (v, u)):       # 原图无向 → 双向
                if a in layers[l] and b in layers[l + 1]:
                    arcs.append(((a, l), (b, l + 1)))
    for l in range(2, L + 1):                   # dst 自环：允许「提前到达后填充」
        arcs.append(((dst, l), (dst, l + 1)))
    return arcs


def maxflow(arcs, s, t):
    """Edmonds-Karp，单位容量"""
    cap = {}
    for a, b in arcs:
        cap[(a, b)] = cap.get((a, b), 0) + 1
    adj = {}
    for a, b in cap:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)

    flow = 0
    while True:
        par, q = {s: None}, deque([s])
        while q and t not in par:
            u = q.popleft()
            for v in adj.get(u, []):
                if v not in par and cap.get((u, v), 0) - cap.get((v, u), 0) > 0:
                    par[v] = u
                    q.append(v)
        if t not in par:
            return flow
        v = t
        while v != s:
            u = par[v]
            cap[(u, v)] -= 1
            cap[(v, u)] = cap.get((v, u), 0) + 1
            v = u
        flow += 1


# ------------------------------------------------ 真实要求（暴力枚举）
def true_disjoint_paths(edges, L, src, dst):
    adj = {}
    for u, v in edges:
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)

    def walk(u, seen, used):
        if len(used) > L:
            return
        if u == dst:
            yield tuple(used)
            return
        for v in adj.get(u, []):
            e = frozenset((u, v))
            if e in used or v in seen:
                continue
            seen.add(v)
            used.append(e)
            yield from walk(v, seen, used)
            used.pop()
            seen.discard(v)

    paths = list(walk(src, {src}, []))
    for k in range(len(paths), 0, -1):
        for comb in itertools.combinations(paths, k):
            es = [e for p in comb for e in p]
            if len(set(es)) == len(es):          # 物理边互不相同
                return k, len(paths)
    return 0, len(paths)


def check(edges, L, src='s', dst='t', tag=''):
    mf = maxflow(build_layered(edges, L, src, dst), (src, 1), (dst, L + 1))
    tr, npath = true_disjoint_paths(edges, L, src, dst)
    ok = mf <= tr
    print(f'{tag} |V|={len(set(sum(([u, v] for u, v in edges), [])))} '
          f'|E|={len(edges)} L={L} -> max-flow={mf}, true={tr}, '
          f'{"一致" if ok else "★ 不一致(算法过判)"}')
    return mf, tr


if __name__ == '__main__':
    print('=' * 68)
    print('例 1：手造反例（s-x-a-b-t 链 + 捷径 s-a）')
    mf, tr = check(DEFAULT_EDGES, 4, tag='  ')
    print(f'  → 算法认为有 {mf} 条不相交路径(判可行)；实际只有 {tr} 条。')
    if mf > tr:
        print('  → 两条"不相交"路径分别在 layer 2 和 layer 3 复用了物理边 a-b / b-t。')

    print('=' * 68)
    print('随机小图统计：看这是不是个例')
    random.seed(0)
    miss = tot = 0
    for _ in range(300):
        n = random.randint(5, 7)
        names = list('sabcdef'[:n])
        names = ['s'] + names[1:-1] + ['t']
        cand = list(itertools.combinations(names, 2))
        E = [e for e in cand if random.random() < 0.45]
        if len(E) < 5 or ('s', 't') in E:
            continue
        L = random.choice([3, 4])
        try:
            mf = maxflow(build_layered(E, L, 's', 't'), ('s', 1), ('t', L + 1))
            tr, _ = true_disjoint_paths(E, L, 's', 't')
        except Exception:
            continue
        tot += 1
        if mf > tr:
            miss += 1
    print(f'  有效样本 {tot} 个，其中 max-flow > 真实值（过判）{miss} 个 '
          f'({100.0*miss/max(tot,1):.1f}%)')
    print('  过判样本的共同点：分层图上的两条流路径在不同层复用了同一条物理边。')
