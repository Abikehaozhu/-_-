import cv2
import numpy as np


def seed_grow(img, output, x, y, thresh):
    """
    对图片进行种子生长
    :param output: 输出的数组，需预先设置好
    :param img:输入的图片
    :param x:输入的x坐标，对应图像的行
    :param y: 输入的y坐标，即x行y列
    :param thresh:种子生长的阈值
    :return:返回种子生长的结果图
    """
    direction = [[0, 1], [1, 1], [1, 0], [1, -1],
                 [0, -1], [-1, -1], [-1, 0], [-1, 1]]  # 八个方向
    stack_x = []  # 存的栈
    stack_y = []
    stack_x.append(x)  # 种子点的坐标进栈
    stack_y.append(y)
    seed_num = int(img[x][y])  # 读取该点的数值
    while len(stack_x) != 0:
        temp_x = stack_x.pop()  # 弹出栈顶元素
        temp_y = stack_y.pop()
        for k in range(8):  # 在8个方向上查找种子点
            current_x = temp_x + direction[k][0]
            current_y = temp_y + direction[k][1]
            if current_x >= img.shape[0] or current_x < 0 \
                    or current_y >= img.shape[1] or current_y < 0:
                continue
            temp_num = int(img[current_x][current_y])  # 读取现在的点的值
            if np.abs(temp_num-seed_num) < thresh and output[current_x][current_y] == 0:
                stack_x.append(current_x)
                stack_y.append(current_y)
                output[current_x][current_y] = 255
    # return np.array(output, dtype=np.uint8)


if __name__ == "__main__":
    img = cv2.imread("img/007.tiff", 0)
    output = np.array([[0] * img.shape[1]] * img.shape[0],dtype=np.uint8)  # 输出的数组
    seed_pos = cv2.imread("result/result.bmp", 0)
    # for x in range(seed_pos.shape[0]):
    #     for y in range(seed_pos.shape[1]):
    #         if seed_pos[x][y] == 255:
    #             print(x, y)
    seed_grow(img, output, 604, 335, 50)
    seed_grow(img, output, 258, 627, 40)
    seed_grow(img, output, 365, 417, 40)
    seed_grow(img, output, 355, 417, 50)
    cv2.imshow("out", output)
    cv2.waitKey(0)
    # cv2.imwrite("result/final_result.bmp", output)
