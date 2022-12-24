import cv2
import time
import numpy as np

pi = np.pi


def circular_filter(img, filter_r, n=40, filter_thresh=0.4):
    """
    对输入的图像进行圆周滤波，查找飞机所在的坐标位置
    :param img: 输入的单通道图像
    :param filter_r: 圆周滤波的半径
    :param n: 圆周率波的点数
    :param filter_thresh: 滤波的阈值，高于阈值则视作存在飞机
    :return: 得到飞机坐标图，二值
    """
    start_time = time.time()
    height = img.shape[0]  # 读取图像的长、宽
    length = img.shape[1]
    output = np.array([[0] * length] * height)  # 输出的矩阵
    for x in range(filter_r, height - filter_r):
        for y in range(filter_r, length - filter_r):  # 逐点进行圆周滤波
            a = 0.0
            b = 0.0
            for i in range(n):  # 进行n点圆周滤波
                a += img[int(x + filter_r * np.sin(2 * pi * i / n)), \
                         int(y + filter_r * np.cos(2 * pi * i / n))] \
                     * np.cos(8 * pi * i / n)
                b += img[int(x + filter_r * np.sin(2 * pi * i / n)), \
                         int(y + filter_r * np.cos(2 * pi * i / n))] \
                     * np.sin(8 * pi * i / n)
            output[x][y] = a ** 2 + b ** 2
    max_num = np.max(output)
    output = output / max_num
    for x in range(height):
        for y in range(length):
            if output[x][y] > filter_thresh:
                output[x][y] = 255
            else:
                output[x][y] = 0
    end_time = time.time()
    print(f"圆周滤波结束，用时{end_time-start_time}s.")
    return output


if __name__ == "__main__":
    img = cv2.imread("img/007.tiff", 0)  # 以灰度图像读取
    plane = circular_filter(img, 20, n=60, filter_thresh=0.4)
    cv2.imwrite("result/result.bmp", plane)
