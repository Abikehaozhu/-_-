#include<stdio.h>
#include<iostream>
#include<math.h>
#include<Windows.h>
#define pi 3.1415926
#define MaxSize 256
using namespace std;
struct HuffmanNode
{
	int weight;
	unsigned char data;
	HuffmanNode* llink;
	HuffmanNode* rlink;
};
struct HuffmanCode
{
	char bit[256]; //���ڻҶ�ͼ�����256��ͳ��ֵ�������256�ı���   
	int len;	// �����ʵ�ʳ���
	unsigned char c;  // ԭʼ�ַ�
};
struct Stack
{
	char Data[MaxSize];
	int index;

	Stack() { this->index = -1; }

	void clear() { this->index = -1; }
	bool isFull() { return this->index >= MaxSize - 1; }
	bool isEmpty() { return this->index == -1; }
	void push(int d) {
		this->Data[++this->index] = d;
	}
	int pop() {
		return this->Data[this->index--];
	}
	int size() { return this->index + 1; }
};
struct MinHeap {

	HuffmanNode* Data[MaxSize];
	int Size;

	bool IsFull()
	{
		return this->Size == MaxSize - 1;
	}

	bool IsEmpty()
	{
		return this->Size == 0;
	}

	bool Push(HuffmanNode* X)
	{
		if (IsFull()) return false;

		int i;
		++this->Size;
		i = this->Size - 1;  //iָ��������е����һ��Ԫ�ص�λ��
		for (; this->Data[(i - 1) / 2]->weight > X->weight && i > 0; i = (i - 1) / 2) //�븸���Ƚ�
		{
			this->Data[i] = this->Data[(i - 1) / 2];  //�����ڸ�������븸��㽻��
		}
		this->Data[i] = X;  //�ҵ��˺��ʵ�λ�ã���X����

		return true;
	}

	void PercDown(int p)
	{
		HuffmanNode* X = this->Data[p]; //ȡ��������ֵ

		int parent, child;
		for (parent = p; parent * 2 + 1 < this->Size; parent = child) { //�����ӽ�����ʱ����
			child = parent * 2 + 1;
			if (child != this->Size - 1 && this->Data[child]->weight > this->Data[child + 1]->weight) { // ���ӽ����ڣ���С�����ӽ��
				child++;
			}
			if (X->weight <= this->Data[child]->weight)
				break; // �ҵ��˺��ʵ�λ��
			else
				this->Data[parent] = this->Data[child]; // �³�X
		}
		this->Data[parent] = X;
	}

	void BuildHeap() {
		for (int i = this->Size / 2 - 1; i >= 0; --i) PercDown(i);
	}

	HuffmanNode* Pop()
	{
		HuffmanNode* maxItem = this->Data[0]; //ȡ��������ŵ����ֵ
		HuffmanNode* X = this->Data[this->Size - 1];
		this->Size--;

		int parent, child;
		for (parent = 0; parent * 2 + 1 < this->Size; parent = child) { //�����ӽ�����ʱ����
			child = parent * 2 + 1;
			if (child < this->Size - 1 && this->Data[child]->weight > this->Data[child + 1]->weight) { // ���ӽ����ڣ��Ҵ������ӽ��
				child++;
			}
			if (X->weight <= this->Data[child]->weight)
				break; // �ҵ��˺��ʵ�λ��
			else
				this->Data[parent] = this->Data[child];  // �³�X   
		}
		this->Data[parent] = X;

		return maxItem;
	}
};
void deleteHuffmanTree(HuffmanNode* T);
HuffmanNode* createHuffmanTree(unsigned char buff[256], int weightArray[256], int n);
//����huffman����weightArray����Ȩ�����飬n�Ǵ���������г���
Stack codeStack;
int huffmanCoding(HuffmanNode* T, HuffmanCode* codeTable, int& index);
void calculateGrayLevelWeight(unsigned char* imgBuf, int width, int height, unsigned char grayLevel[],
	int grayLevelWeight[], int* grayLevelCount);
void generateFile(unsigned char* imgBuf, int width, int height, const char* fileout);//���������벢д���ļ�
int decode8BitBuf(char codebuf, int ncode, HuffmanCode* codeTable, int matchPos[],
	unsigned char decodedBuf[], int decodedNum, int imgSize);
unsigned char* decoding(const char* filename, int* width, int* height);
bool saveBmp(const char* bmpName, unsigned char* imgBuf, int width, int height, int byteCount);
struct point//�����
{
	int x;
	int y;
};
void CircularFilter(int input[], double output[], int m, int n, 
	int filter_r, int N, double filter_thresh);//�ԻҶ�ͼ�����Բ���˲�
int find_seed(double filter_out[], int m, int n, point p[],
	int  cluster_r, int  cluster_r_y);//Ѱ�����ӵ�
void regionGrowing_raw(unsigned char* buf, int m,
	int n, int seed_y, int seed_x, int thresh, unsigned char* bufOut);//����������
unsigned char* readBmp(const char* bmpName, int* width, int* height, int* byteCount);//����bmp�ļ�




int main()
{
	//�������趨������007
	int m = 985, n = 769;// ͼƬ�ĳ���Ϊm=985,n=769
	int filter_r = 20, N = 60;//Բ���˲��İ뾶�Լ�����
	double filter_thresh = 0.4;//Բ���˲�����ֵ
	int cluster_r_x = 4, cluster_r_y = 4;//ѡȡ���ӵ�Ԥ�������x��y���
	int thresh = 70;//������������ֵ


	//�����ļ��Ķ�ȡ
	FILE* fp;
	errno_t err = fopen_s(&fp, "data\\gray.raw", "r");//�ļ������ô���Ϊ�Ҷ�ͼ����PS�����raw��ʽ
	if (err != 0)
	{//����ļ��Ƿ��
		printf("siu_fail");
		return 0;
	}
	unsigned char* bufGray = new unsigned char[m * n];//���ڱ���Ҷ�ͼ��
	fread(bufGray, sizeof(unsigned char), m * n, fp);//��ͼ�����bufGray��
	fclose(fp);
	int* buf = new int[m * n];//���޷����ַ���ת��Ϊint���͵ľ���
	double* filter_result = new double [m*n];//Բ���˲���ľ���
	for (int i = 0; i < m * n; i++)
	{//���о������͵�ת���Ѿ�filter_result�ĳ�ʼ��
		buf[i] = (int)bufGray[i];
		filter_result[i] = 0.0;
	}
	CircularFilter(buf, filter_result,m,n,filter_r,N,filter_thresh);//����Բ���˲�
	printf("Բ���˲����\n");
	point* position = new point[10000];//���������������
	int seed_num=find_seed(filter_result, m, n, position, cluster_r_x,cluster_r_y);//Ѱ������������
	unsigned char* bufOut = new unsigned char[m * n];//���������
	for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				bufOut[i * m + j] = 255;
			}
		}
	for (int i = 0; i < seed_num; i ++)
	{//��Ѱ�ҵ������ӵ������������
		regionGrowing_raw(bufGray, m, n, position[i].x+3, 
			position[i].y-18, thresh, bufOut);
	}
	printf("�����������\n");
	//err = fopen_s(&fp, "output.raw", "wb");
	//fwrite(bufOut, sizeof(unsigned char), m * n, fp);
	//fclose(fp);
	char fn[255] = "data\\outcome.huff";
	generateFile(bufOut, m, n, fn);//�����յĽ��д��huff�ļ�������
	printf("�Ѿ�д��huff�ļ�\n");

	/*char fy[255] = "data\\find_007.bmp";
	saveBmp(fy, bufOut, m, n, 1);
	printf("�Ѿ�����Ϊbmp�ļ�\n");*/
	/*�������huff�ļ��Ľ���д��bmp�ļ�*/
	printf("�Ƿ��huff�ļ����н��룿(Y/N):\n");
	char siusiu;
	cin >> siusiu;
	if (siusiu == 'Y')
	{//�����Ҫ���н���д��bmp
		char fx[255] = "data\\recover.bmp";
		int width, height, byteCount = 1;
		//����
		unsigned char* buf2 = decoding(fn, &width, &height);
		printf("�������\n");
		saveBmp(fx, buf2, width, height, byteCount);
		printf("�Ѿ�дΪbmp\n");
	}

	//���潫�˲����д��txt�ļ���
	
	//err = fopen_s(&fp, "012.txt", "w");//���潫�˲����д��txt�ļ���
	//for (int i = 0; i < n; i++)
	//{
	//	for (int j = 0; j < m; j++)
	//	{
	//		fprintf(fp, "%d%c", int(bufOut[m * i + j]), ' ');
	//	}
	//	fprintf(fp, "%c",'\n');
	//}
	//fclose(fp);
	
	delete[]bufGray;
	delete[]buf;
	delete[]filter_result;
	delete[]position;
	delete[]bufOut;
	
}




void CircularFilter(int input[], double output[], int m,
	int n,int filter_r, int N, double filter_thresh)
{	//��������Ϊ�����ͼ�������ͼ��ͼ�񳤡�
	//ͼ���Բ���˲��뾶��ѡȡ�ĵ������˲���ֵ
	long double A, B;
	double max = -1;//��ȡ���ֵ�������һ��
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			A = 0,B = 0;
			if (i - filter_r >= 0 && i + filter_r < n &&
				j - filter_r >= 0 && j + filter_r < m)
			{//�ж��Ƿ���ͼ����
				for (int k = 0; k < N; k++)
				{//����Բ���˲�
					A = A + input[m*(int(i + filter_r * sin(2 * pi * k / N)))+
						int(j + filter_r * cos(2 * pi * k / N))] * cos(8 * pi * k / N);
					B = B + input[m*(int(i + filter_r * sin(2 * pi * k / N)))+
						int(j + filter_r * cos(2 * pi * k / N))] * sin(8 * pi * k / N);
				}
				output[i * m + j] = pow(A,2) + pow(B, 2);
				if (output[i * m + j] > max)//��ͼ����ȡ���ֵ
					max = output[i * m+ j];
			}
		}
	}
	for (int i = 0; i < m * n; i++)
	{//���˲����ͼ����ж�ֵ��
		output[i] = output[i] / max;//��ͼ����й�һ��
		if (output[i] < filter_thresh)//С��Ԥ������ֵʱΪ��
			output[i] = 0;
		else
			output[i] = 255;
	}
}
int find_seed(double filter_out[], int m,int n,point p[], 
	int  cluster_r_x,int  cluster_r_y)
{	//Ѱ�����ӵ�,��������Ϊ������˲������
	//���󳤡��������x�������y������󡢾���İ뾶
	int top = -1;
	for (int i = cluster_r_x; i < n-cluster_r_x; i++)
	{
		for (int j = cluster_r_y; j < m-cluster_r_y; j++)
		{
			if (filter_out[(i - cluster_r_x)* m + j] != 0
				&& filter_out[(i + cluster_r_x) * m + j] != 0
				&& filter_out[i * m  + j-cluster_r_y] != 0 
				&& filter_out[i * m + j + cluster_r_y] != 0)
			{
				top++;
				p[top].x = i;
				p[top].y = j;
			}
		}
	}
	return top + 1;//�������ӵ�����Ĵ�С
}
void regionGrowing_raw(unsigned char* buf, int m, int n, 
	int seed_y, int seed_x, int thresh, unsigned char* bufOut)
{//�������������������޷����ַ��������������,��������Ϊ�����ͼ��
 //ͼ�񳤡�ͼ������ӵ��xy���ꡢ������ֵ���������
	int i, j, k;//ѭ������
	//��ά����direction�������ĵ�8����������õ���x��y�����ϵ�ƫ��,
	//���е�һάΪx�����ƫ��,�ڶ�άΪy�����ƫ��
	int direction[8][2] = { { 0,1 },{ 1,1 },{ 1,0 },{ 1,-1 },
		{ 0,-1 },{ -1,-1 },{ -1,0 },{ -1,1 } };

	//ջ���룬������ջ��ָ��
	point* stack = new point[m * n];
	int top = -1;//ջ��ָ��

	//��ǰ������ĵ�͵����ĵ�
	int currentPoint_x, currentPoint_y, popPoint_x, popPoint_y;

	int temp1, temp2;//��ʱ����,������ӵ��Լ����������ֵ

	//��¼���ӵ��ֵ
	temp1 = *(buf + seed_y * m + seed_x);

	//���������ӵ��ñ��0,��ջ
	bufOut[seed_y * m + seed_x] = 0;
	top++;
	stack[top].x = seed_x;
	stack[top].y = seed_y;

	//��ջ
	while (top > -1) {
		//����ջ��Ԫ��,��Ԫ���Ѿ�������
		popPoint_x = stack[top].x;
		popPoint_y = stack[top].y;
		top--;

		//���챻����Ԫ����Χ�Ƿ���û��������Ԫ��
		for (k = 0; k < 8; k++) {
			currentPoint_x = popPoint_x + direction[k][0];
			currentPoint_y = popPoint_y + direction[k][1];
			if (currentPoint_x >= m || currentPoint_x < 0 || currentPoint_y >= n || currentPoint_y < 0)
			{
				continue;
			}
			temp2 = buf[currentPoint_y * m + currentPoint_x];
			if (abs(temp1 - temp2) < thresh && bufOut[currentPoint_y * m + currentPoint_x] == 255)
			{
				top++;
				stack[top].x = currentPoint_x;
				stack[top].y = currentPoint_y;
				bufOut[currentPoint_y * m + currentPoint_x] = 0;
			}

		}//k

	}//while(top)

	 //���������
	delete[]stack;
}
void deleteHuffmanTree(HuffmanNode* T)
{
	if (T == 0)//����������
		return;
	struct HuffmanNode* s[200], * t;
	int top;
	top = 0;
	s[top] = T;//������ջ
	while (top > -1) {
		t = s[top];
		top--;
		if (t->rlink != 0) {//��������ջ
			top++;
			s[top] = t->rlink;
		}
		if (t->llink != 0) {//��������ջ
			top++;
			s[top] = t->llink;
		}
		delete t;//�ͷŽ��
	}
}
HuffmanNode* createHuffmanTree(unsigned char buff[256], int weightArray[256], int n)
{
	MinHeap heap;
	heap.Size = n;
	for (int i = 0; i < n; ++i)
	{
		heap.Data[i] = new HuffmanNode;
		heap.Data[i]->data = buff[i];
		heap.Data[i]->llink = nullptr;
		heap.Data[i]->rlink = nullptr;
		heap.Data[i]->weight = weightArray[i];
	}
	heap.BuildHeap();

	for (int i = 0; i < n - 1; ++i)
	{
		HuffmanNode* node = new HuffmanNode;
		node->llink = heap.Pop();
		node->rlink = heap.Pop();
		node->weight = node->llink->weight + node->rlink->weight;
		heap.Push(node);
	}
	return heap.Pop();
}
bool saveBmp(const char* bmpName, unsigned char* imgBuf, int width, int height, int byteCount)
{
	if (!imgBuf)
		return 0;

	//�Ҷ�ͼ����ɫ��ռ�1024����ɫͼ��û����ɫ��
	int colorTable = 0;
	if (byteCount == 1) colorTable = 1024;

	//һ�������ֽ���Ϊ4�ı���
	int lineByte = (width * byteCount + 3) / 4 * 4;

	FILE* fp;
	fopen_s(&fp, bmpName, "wb");
	if (fp == 0) return 0;

	//��д�ļ�ͷ
	BITMAPFILEHEADER fileHead;
	fileHead.bfType = 0x4D42;
	fileHead.bfSize =
		sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + colorTable + lineByte * height;
	fileHead.bfReserved1 = 0;
	fileHead.bfReserved2 = 0;
	fileHead.bfOffBits = 54 + colorTable;
	fwrite(&fileHead, sizeof(BITMAPFILEHEADER), 1, fp);

	// ��д��Ϣͷ
	BITMAPINFOHEADER head;
	head.biBitCount = byteCount * 8;
	head.biClrImportant = 0;
	head.biClrUsed = 0;
	head.biCompression = 0;
	head.biHeight = height;
	head.biPlanes = 1;
	head.biSize = 40;
	head.biSizeImage = lineByte * height;
	head.biWidth = width;
	head.biXPelsPerMeter = 0;
	head.biYPelsPerMeter = 0;
	fwrite(&head, sizeof(BITMAPINFOHEADER), 1, fp);

	//��ɫ����  
	if (colorTable == 1024)
	{
		unsigned char table[1024];
		for (int i = 0; i < 256; i++)
		{
			*(table + i * 4 + 0) = i;
			*(table + i * 4 + 1) = i;
			*(table + i * 4 + 2) = i;
			*(table + i * 4 + 3) = 0;
		}
		fwrite(table, 1024, 1, fp);
	}

	//׼�����ݲ�д�ļ�
	unsigned char* buf = new unsigned char[height * lineByte];
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * byteCount; j++)
			*(buf + i * lineByte + j) = *(imgBuf + i * width * byteCount + j);
	}
	fwrite(buf, height * lineByte, 1, fp);

	delete[]buf;

	fclose(fp);

	return 1;
}
int decode8BitBuf(char codebuf, int ncode, HuffmanCode* codeTable, int matchPos[],
	unsigned char decodedBuf[], int decodedNum, int imgSize)
{
	int i, j;
	int t, k, l;//
	for (i = 0; i < 8; i++)   //���λ��ÿ������
	{
		t = (codebuf & 1 << i) != 0;//��õ�ǰ����

		for (j = 0; j < ncode; j++) {//����е�j���ַ�������
			if (matchPos[j] == -1)//�Ѿ�ʧ�ܣ��������
				continue;

			if (codeTable[j].bit[matchPos[j]] == t) {//����д�ƥ���λ���뵱ǰ������ͬ				
				matchPos[j]++;
			}
			else {//ûƥ����,��-1��������һ��
				matchPos[j] = -1;
				continue;
			}

			if (matchPos[j] == codeTable[j].len) {//����һ���ַ�,һ�������볤���������������ĳ���
				decodedBuf[decodedNum] = codeTable[j].c;
				decodedNum++;

				//matchPosҪ���³�ʼ��������������ƥ��
				for (l = 0; l < ncode; l++)
					matchPos[l] = 0;

				break;//��ƥ����һ���ַ�ʱ������������ַ��Ͳ���ƥ�䣬����ѭ��
			}
		}//j

		if (decodedNum == imgSize)
			break;
	}//i
	return decodedNum;//��ǰ�ѽ������
}
unsigned char* decoding(const char* filename, int* width, int* height)
{
	// ��ҵ
	FILE* fp;
	fopen_s(&fp, filename, "rb");
	//fp = fopen(filename, "rb");
	int ncode;//�����
	fread(&ncode, sizeof(int), 1, fp);
	HuffmanCode huffmanCode[256];
	for (int i = 0; i < ncode; i++)//����Huffman��
	{
		fread(&huffmanCode[i].c, sizeof(huffmanCode[i].c), 1, fp);
		fread(&huffmanCode[i].len, sizeof(huffmanCode[i].len), 1, fp);
		fread(&huffmanCode[i].bit, sizeof(char), huffmanCode[i].len, fp);
	}
	fread(width, sizeof(int), 1, fp);
	fread(height, sizeof(int), 1, fp);
	int w = *width;
	int h = *height;
	unsigned char* img_buf = new unsigned char[h * w];
	int matchPos[256] = { 0 };
	int decodedNum = 0;
	while (decodedNum < h * w)//���ν���ÿһ��Huffman���벢����buf��
		//ÿ�ζ�8λ�������뺯��
	{
		char readbit;
		fread(&readbit, sizeof(char), 1, fp);
		//			printf("readbit: %d ",readbit);
		decodedNum = decode8BitBuf(readbit, ncode, huffmanCode, matchPos, img_buf, decodedNum, (*width) * (*height));
		//			printf("decodenum:%d\n",decodedNum);
	}

	//	return nullptr;
	return img_buf;
}
void generateFile(unsigned char* imgBuf, int width, int height, const char* fileout)
{
	//ͳ��ͼ����Ƶ�ʷ�0�ĻҶȼ�����Ȩ��
	unsigned char grayLevel[256];
	int grayLevelWeight[256], grayLevelCount;
	calculateGrayLevelWeight(imgBuf, width, height, grayLevel, grayLevelWeight, &grayLevelCount);

	//�����������������ɱ�������
	HuffmanNode* HuffTree = createHuffmanTree(grayLevel, grayLevelWeight, grayLevelCount);//����huffman��
	HuffmanCode huffmanCode[256];//��Ź���������
	codeStack.clear();
	int index = 0;
	int ncode = huffmanCoding(HuffTree, huffmanCode, index);//����n��Ҷ��㣬huffmanTree����һ����2n-1����

	//�����ļ�
	FILE* fp;
	fopen_s(&fp, fileout, "wb");

	fwrite(&ncode, sizeof(ncode), 1, fp);
	for (int i = 0; i < ncode; ++i)
	{
		fwrite(&huffmanCode[i].c, sizeof(huffmanCode[i].c), 1, fp);
		fwrite(&huffmanCode[i].len, sizeof(huffmanCode[i].len), 1, fp);
		fwrite(&huffmanCode[i].bit, sizeof(char), huffmanCode[i].len, fp);
	}

	//ͼ�����д���ļ�������ʱ����������ʵ����ݿռ�
	fwrite(&width, sizeof(int), 1, fp);
	fwrite(&height, sizeof(int), 1, fp);

	//��ÿһ���ַ�����д��Ĺ���
	int i, j, k, t1, t2;
	char codebuf;//���뻺����
	t1 = 0;//�������ڱ���λ��
	for (i = height - 1; i >= 0; i--) {
		for (j = 0; j < width; j++) {
			//�ҵ���Ӧ�Ҷȼ���Ȩ��
			for (k = 0; k < ncode; k++) {
				if (huffmanCode[k].c == imgBuf[i * width + j])
					break;
			}
			for (t2 = 0; t2 < huffmanCode[k].len; t2++) {//������������д�뻺����
				if (huffmanCode[k].bit[t2] == 0)
					codebuf &= ~(1 << t1);//����λ����0
				else
					codebuf |= (1 << t1);//����λ����1
				t1++;
				if (t1 == 8) {//���뻺��������д���ļ�
					fwrite(&codebuf, sizeof(char), 1, fp);
					t1 = 0;
				}
			}
			if (i == height - 1 && j == width - 1) {//�������һ������
				if (t1 != 0)//���һ�λ�������û������д���ļ�������ѭ��
					fwrite(&codebuf, sizeof(char), 1, fp);
			}
		}//j
	}//i
	fclose(fp);

	deleteHuffmanTree(HuffTree);
}
int huffmanCoding(HuffmanNode* T, HuffmanCode* codeTable, int& index)
{
	if (T->llink && T->rlink)  //Huffman�������ڶ�Ϊ1�Ľ��
	{
		codeStack.push(0);
		huffmanCoding(T->llink, codeTable, index);
		codeStack.push(1);
		huffmanCoding(T->rlink, codeTable, index);
	}
	else
	{
		codeTable[index].len = codeStack.size();
		codeTable[index].c = T->data;
		for (int i = 0; i <= codeStack.index; ++i)
			codeTable[index].bit[i] = codeStack.Data[i];
		index++;
	}

	codeStack.pop();

	return index;	//����code����Ŀ
}
void calculateGrayLevelWeight(unsigned char* imgBuf, int width, int height, unsigned char grayLevel[],
	int grayLevelWeight[], int* grayLevelCount)
{
	int i, j, hist[256] = { 0 };
	for (i = 0; i < height; i++) {//ͳ��ֱ��ͼ
		for (j = 0; j < width; j++) {
			hist[*(imgBuf + i * width + j)]++;
		}
	}
	for (i = 0; i < 256; i++)//ͳ��Ƶ��/Ȩ��
		grayLevelWeight[i] = 0;
	int count = 0;
	for (i = 0; i < 256; i++) {//ͳ��Ƶ�ʲ�Ϊ0�ĻҶȼ����������ɹ�������
		if (hist[i] != 0) {
			grayLevel[count] = i;
			grayLevelWeight[count] = hist[i];
			count++;
		}
	}
	*grayLevelCount = count;//ʵ�ʻҶȼ���Ŀ
}
unsigned char* readBmp(const char* bmpName, int* width, int* height, int* byteCount)
{
	FILE* fp;
	fopen_s(&fp, bmpName, "rb");
	if (fp == 0) return 0;
	fseek(fp, sizeof(BITMAPFILEHEADER), 0);

	int w, h, b;
	BITMAPINFOHEADER head;
	fread(&head, sizeof(BITMAPINFOHEADER), 1, fp);
	w = head.biWidth;
	h = head.biHeight;
	b = head.biBitCount / 8;
	int lineByte = (w * b + 3) / 4 * 4;

	if (b == 1)
		fseek(fp, 1024, 1);
	unsigned char* imgBuf = new unsigned char[w * h * b];
	for (int i = 0; i < h; i++)
	{
		fread(imgBuf + i * w * b, w * b, 1, fp);
		fseek(fp, lineByte - w * b, 1);
	}
	fclose(fp);

	*width = w, * height = h, * byteCount = b;

	return imgBuf;
}