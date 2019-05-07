from __future__ import absolute_import, division, print_function
import numpy as np
import random
import pywt
import cv2
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import util


# supported transforms
TRANSFORM_MODES = ['fwht', 'fft']


# check if x is a power of 2
def is_binary(x):
    lvl = 0
    while (1 << lvl) < x: lvl += 1
    return (1 << lvl) == x


# Fast Hahamard Transform without normalization
def fht(X):
    Y = np.array(X, copy=True)
    N = X.shape[0]
    lvl = 0
    while (1 << lvl) < N:
        mask = np.zeros((N), dtype=bool)
        mask[np.arange(N) & (1 << lvl) == 0] = True
        Yl = Y[mask]
        Yr = Y[~mask]
        Y[mask] = Yl + Yr
        Y[~mask] = Yl - Yr
        lvl += 1

    return Y


# Fast Walsh-Hadamard Transform
# from https://github.com/dingluo/fwht/blob/master/FWHT.py
def fwht(X, normalize=True):
    N = X.shape[0]
    Nrest = X.shape[1:]
    G = N//2 # Number of Groups
    M = 2 # Number of Members in Each Group

    # step first leevl
    X_next = np.zeros((M, G) + Nrest)
    X_next[0, :] = X[0::2] + X[1::2]
    X_next[1, :] = X[0::2] - X[1::2]
    X = X_next

    while X.shape[1] > 1:
        # extract blocks
        even_l = X[0:M:2, 0:G:2]
        even_r = X[0:M:2, 1:G:2]
        odd_l = X[1:M:2, 0:G:2]
        odd_r = X[1:M:2, 1:G:2]
        G //= 2
        M *= 2

        # step another level
        X_next = np.zeros((M, G) + Nrest)
        X_next[0:M:4, 0:G] = even_l + even_r
        X_next[1:M:4, 0:G] = even_l - even_r
        X_next[2:M:4, 0:G] = odd_l - odd_r
        X_next[3:M:4, 0:G] = odd_l + odd_r
        X = X_next

    X = X[:, 0]
    if normalize: X /= np.sqrt(N)
    return X


# 2D Fast Walsh-Hahamard Transform without normalization
def fwht2(X):
    return fwht(fwht(X).T).T


# dft + shift
def dft_shifted(X):
    dft = cv2.dft(X.astype(np.float32), flags=cv2.DFT_COMPLEX_OUTPUT)
    dft = np.fft.fftshift(dft)
    return dft


"""
Feature extraction by multilevel sampling
    transform   'fwht': Fast Walsh-Hadamard Transform
                'fft':  Fast Fourier Transform
    seed        random generator seed number
"""
class MultilevelRandomFeature:

    def __init__(self, img_shape, s_ratio, transform_mode='fwht', seed=None):
        # argument validatation
        assert 0.0 <= s_ratio <= 1.0
        assert transform_mode in TRANSFORM_MODES
        if transform_mode in ['fwht', 'fft']:
            assert is_binary(img_shape[0]) and is_binary(img_shape[1])
            assert img_shape[0] == img_shape[1]

        # save state
        self.img_shape = img_shape
        self.s_ratio = s_ratio
        self.n_feature = int(np.prod(img_shape) * s_ratio)
        self.transform_mode = transform_mode
        self.seed = seed
        self.rng = random.Random(self.seed)
        self.rng_np = np.random.RandomState(self.seed)

        # generate multilevel sampling indices
        self.init_indices()


    # initialize indices accordingly
    def init_indices(self):
        if self.n_feature == np.prod(self.img_shape):
            self.indices = np.meshgrid(np.arange(self.img_shape[0]),
                                       np.arange(self.img_shape[1]))
        elif self.transform_mode == 'fwht':
            scheme = self.uniform_rect_samp_scheme()
            pos = self.random_rect_subsamp(scheme)  # 2 x M indices
            self.indices = (pos[0], pos[1])
        elif self.transform_mode == 'fft':
            nu = (np.log2(min(self.img_shape))).astype(int)
            bounds = [2**(nu-4), 2**(nu-2), 2**(nu-1)]
            pos = self.random_circ_subsamp(bounds)  # 2 x M indices
            self.indices = (pos[0], pos[1])


    # extract feature from img
    def transform(self, img):
        img = self.fit_image(img)
        if self.transform_mode == 'fwht':
            Y = fwht2(img)  # transform

        elif self.transform_mode == 'fft':
            Y = dft_shifted(img)  # transform
        Y = Y[self.indices]  # sample
        return Y


    # fit image to valid size
    def fit_image(self, img):
        new_img = np.zeros(self.img_shape)
        if img.shape != self.img_shape:
            new_img = cv2.resize(img, dsize=self.img_shape)
        return new_img


    # generate rectangle sampling scheme
    def uniform_rect_samp_scheme(self):
        nu = int(np.log2(self.img_shape[0]))
        scheme = np.ones((nu,), dtype=int)

        # scale by 2 each level
        i = 0
        while np.sum(scheme**2) <= self.n_feature/3:
            scheme[i+1:nu] += 2**i
            i += 1
        scheme[i:nu] = 2**(i-1)
        
        # fill up
        while np.sum(3 * (scheme**2)) <= self.n_feature:
            scheme[i:nu] += 1
        scheme[i:nu] -= 1
        scheme = 3 * (scheme**2)
        scheme[0] = 4
        
        # add residual in round robin fashion
        k = 0
        while np.sum(scheme) < self.n_feature:
             s = i + (k % (nu-i))
             scheme[s] += 1
             k += 1

        return scheme


    # sample with rectangle sampling scheme
    def random_rect_subsamp(self, scheme):
        b = 2**(np.arange(len(scheme)+1))
        pos_buffer = []
        
        if scheme[0] != 4:
            raise ValueError('First scheme element must be 4')
        pos_buffer.append([0, 0])
        pos_buffer.append([0, 1])
        pos_buffer.append([1, 0])
        pos_buffer.append([1, 1])
    
        # sample all position
        k = 1
        while b[k+1]**2 - b[k]**2 == scheme[k]:
            low = b[k]
            high = b[k+1]
            for i in range(low, high):
                for j in range(0, high):
                    pos_buffer.append([i, j])

                for j in range(0, low):
                    pos_buffer.append([j, i])
            k += 1
        
        # random sample the rest
        while k < len(scheme):
            low = b[k]
            high = b[k+1]
            for _ in range(scheme[k]):
                elem = ((high-1) * self.rng_np.rand(2) + 1).astype(int)
                while np.max(elem) < low:
                    elem = ((high-1)*self.rng_np.rand(2) + 1).astype(int)
                pos_buffer.append(elem)
            k += 1
        return np.array(pos_buffer).T


    # sample multiple circular levels
    def random_circ_subsamp(self, bounds):
        raise NotImplementedError('work in progess...')


if __name__ == '__main__':
    # demo for MultilevelRandomFeatures
    mlrf = MultilevelRandomFeature(
        img_shape= (256, 256),  # image size (fwht requires 2**n)
        s_ratio= 0.2,  # sparsity ratio
        transform_mode= 'fwht',  # transform function
        seed= 42,  # random state
    )

    img_orig = cv2.imread('/path/to/image/pic.png', 0)  # YOUR IMAGE HERE
    img = mlrf.fit_image(img_orig)
    img = (img - img.min()) / (img.max() - img.min())
    Y = fwht2(img)
    img2 = fwht2(Y)

    # plot focus on low frequency contents
    NCUT = 64
    fig, axes = plt.subplots(ncols=2)
    axes[0].imshow(Y[:, :])
    axes[1].imshow(Y[:NCUT, :NCUT])
    rec = patches.Rectangle((0, 0), NCUT, NCUT, linewidth=0.5, 
                            edgecolor='r', facecolor='none')
    axes[0].add_patch(rec)

    # test reconstruction from xs
    fig, axes = plt.subplots(ncols=3)
    Y2 = np.zeros(Y.shape)
    Y2[tuple(mlrf.indices)] = Y[tuple(mlrf.indices)]
    axes[0].imshow(img)
    axes[1].imshow(Y)
    axes[2].imshow(fwht2(Y2))

    # plot multiple level details
    nu = int(np.log2(img.shape[0]))
    NC, NR, S = nu+1, 2, 1.5
    fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(S*NC, S*NR))
    for ax1, ax2, p in zip(axes[0], axes[1], range(nu+1)):
        Y3 = np.array(Y2)
        Y3[2**p:, :] = 0.0
        Y3[:, 2**p:] = 0.0
        ax1.imshow(fwht2(Y3))  # show accumulated image
        ax1.axis('off')
        if p > 0:
            Y3[:2**(p-1), :] = 0.0
            Y3[:, :2**(p-1)] = 0.0
        ax2.imshow(fwht2(Y3))  # show additional content
        ax2.axis('off')
    plt.show()

    # plot of rmses per sparsity
    def mlrf_sparsity(img, Y, img_shape, s_ratio):
        mlrf = MultilevelRandomFeature(
            img_shape= img_shape,
            s_ratio= s_ratio,
            transform_mode= 'fwht',
            seed= 42)
        Y2 = np.zeros(Y.shape)
        Y2[mlrf.indices] = Y[mlrf.indices]
        nu = int(np.log2(img.shape[0]))
        err_lvl = []
        for p in range(nu):
            Y3 = np.array(Y2)
            Y3[2**p:, :] = 0.0
            Y3[:, 2**p:] = 0.0
            img3 = fwht2(Y3)
            err_lvl.append(((img - img3)**2).sum()**0.5)
        return err_lvl
    rmses = []
    Xs = np.linspace(0, 5, 20)
    Xs = np.exp(Xs) / (np.exp(Xs).max() + 1)
    for s_ratio in Xs:
        print('reconstructing with s_ratio= %f'%(s_ratio))
        rmses.append(mlrf_sparsity(img, Y, (256, 256), s_ratio))
    rmses = np.array(rmses)

    fig, axes = plt.subplots(ncols=2)
    for rmse in rmses:
        axes[0].plot(np.arange(rmse.shape[0]), rmse, 'x-')
    for p, rmse in zip(range(rmses.shape[1]), rmses.T):
        axes[1].plot(Xs, rmse, '-', label=str(2**p))
        axes[1].set_xlabel('sparsity')
        axes[1].set_ylabel('RMSE')
        axes[1].legend()
    plt.show()