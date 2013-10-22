#include "FEL_decoder_gpu.h"

__device__
void flexibleFractalDecoding( float* original, float* decoded, int nBin, int flip, int shift)
{
	int m = 0;
	if (flip == 0)
	{
		for (int i = 0; i < nBin; i++)
		{
			m = i + shift;
			if (m >= nBin) 
				m = m - nBin; 
		
			decoded[m] = original[i];
		}// end for
	}
	else 
	{
		float temp[64];
		for (int i = 0; i < nBin; i++) 
		{
			temp[i] = original[nBin - 1 - i];
		}
		
		for (int i = 0; i < nBin; i++) 
		{
			m = i + shift;
			if (m >= nBin)
				m = m - nBin;
			decoded[m] = temp[i];
		}		
	}
}

__global__ void
d_divideBlock( int divPar)
{
	// TODO: WARNING! here we could not use x in this function, otherwise cuda kernel will have error, the problem occurs when we use span[i].low.y
	//int divPar = 30;
	// NOTE: here we assume the input volume data is nDimension*nDimension*nDimension, which means, we didin't consider dataset such as isabel 500x500x100
	// divPar is the size of each block specified by user
	//printf("start d_divideBlock!\n");
	if (divPar > nDimension) {
		printf("wrong divPar! %d\n", divPar);
	}

	int n = 0; // number of blocks divided in each dimension, so that n*n*n is the total number of blocks	
	int2 spanX[flexNBin];	// for each dimension, record the spanLow and spanHigh, x is spanLow, y is spanHigh
	int2 spanY[flexNBin];
	int2 spanZ[flexNBin];
	if (nDimension % divPar != 0) {
		n = nDimension / divPar + 1;
		for (int i = 0; i < n; i++) {
			spanX[i].x = 1 + i * divPar;
			spanY[i].x = 1 + i * divPar;
			spanZ[i].x = 1 + i * divPar;
			if (i != n - 1) {
				spanX[i].y = (i + 1) * divPar;
				spanY[i].y = (i + 1) * divPar;
				spanZ[i].y = (i + 1) * divPar;
			}
			else {
				spanX[i].y = flexNBin;
				spanY[i].y = flexNBin;
				spanZ[i].y = flexNBin;
			}
		}
	}
	else {
		n = nDimension / divPar;
		for (int i = 0; i < n; i++) {
			spanX[i].x = 1 + i * divPar;
			spanY[i].x = 1 + i * divPar;
			spanZ[i].x = 1 + i * divPar;
			spanX[i].y = (i + 1) * divPar;
			spanY[i].y = (i + 1) * divPar;
			spanZ[i].y = (i + 1) * divPar;
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				flexBlock[i*n*n + j*n + k].low.x = spanX[i].x;
				flexBlock[i*n*n + j*n + k].low.y = spanY[j].x;
				flexBlock[i*n*n + j*n + k].low.z = spanZ[k].x;
				flexBlock[i*n*n + j*n + k].high.x = spanX[i].y;
				flexBlock[i*n*n + j*n + k].high.y = spanY[j].y;
				flexBlock[i*n*n + j*n + k].high.z = spanZ[k].y;
				//printf("blockLow(%2d, %2d, %2d), blockHigh(%2d, %2d, %2d)\n", flexBlock[i*n*n + j*n + k].low.x, flexBlock[i*n*n + j*n + k].low.y, flexBlock[i*n*n + j*n + k].low.z, flexBlock[i*n*n + j*n + k].high.x, flexBlock[i*n*n + j*n + k].high.y, flexBlock[i*n*n + j*n + k].high.z);
			}
		}
	}
	nFlexBlock = n * n * n;
	nFlexBlockX = n;
	nFlexBlockY = n;
	nFlexBlockZ = n;
	printf("nFlexBlock: %d\n", nFlexBlock);
	// TODO: call d_decompose to decompose each span and query the histogram
	//printf("finsh d_divideBlock!\n");
}

__global__ void
d_sumSanHistogram(int n) {
	int cId = blockIdx.x;

	// initialize cornerSumHistogram
	for (int i = 0; i < flexNBin; i++) {
		cornerSumHistogram[cId][i] = 0;
	}

	int totalWeight = 0;
	for (int i = 0; i < nSubSpan[cId]; i++) {
		Span currentSpan = subSpan[cId][i];
		int weight = (currentSpan.high.x - currentSpan.low.x + 1) * (currentSpan.high.y - currentSpan.low.y + 1) * (currentSpan.high.z - currentSpan.low.z + 1);
		for (int j = 0; j < flexNBin; j++) {
			cornerSumHistogram[cId][j] += cornerHistogram[cId][i][j] * weight;
			//if (n == 5 && cornerHistogram[cId][i][j] > 0.000001) {
			//	printf("nSubSpan[%2d], binId[%2d], freq: %f\n", i, j, cornerHistogram[cId][i][j]);
			//}
			//if (n == 5) {
			//	printf("cornerSumHistogram[%2d][%2d]: %f\n", cId, j, cornerSumHistogram[cId][j]);
			//}
		}
		totalWeight += weight;
	}

	/* old normalize method
	// normalize cornerSumHistogram
	float total = 0;
	for (int i = 0; i < flexNBin; i++) {
		total += cornerSumHistogram[cId][i];
	}
	printf("total[%2d]: %f\n", cId, total);
	for (int i = 0; i < flexNBin; i++) {
		cornerSumHistogram[cId][i] = cornerSumHistogram[cId][i] / total;
		if (cornerSumHistogram[cId][i] < 0 || cornerSumHistogram[cId][i] > 1) {
			printf("wrong cornerSumHistogram[%d][%d]: %f\n", cId, i, cornerSumHistogram[cId][i]);
		}
	}*/

	//if (n == 13) {
	//	printf("\nhistogram of block 13:\n");
	//}
	/*
	float total = 0;
	for (int i = 0; i < flexNBin; i++) {
		cornerSumHistogram[cId][i] = cornerSumHistogram[cId][i] / totalWeight;
		//if (n == 13) {
		//	printf("%1.6f\t", cornerSumHistogram[cId][i]);
		//}
		total += cornerSumHistogram[cId][i];
	}
	if (total < 0.999999 || total > 1.000001) {
		printf("n = %2d, corner = %d, total = %f\n", n, cId, total);
	}*/
}

__global__ void
d_computeBlock(int n) {
	// n is the index of current block
	float blockHistogram[flexNBin];
	
	//if (n == 13) {
	//	Span cs = flexBlock[n];
	//	printf("histogram of block 13: spanLow(%d, %d, %d), spanHigh(%d, %d, %d)\n", cs.low.x, cs.low.y, cs.low.z, cs.high.x, cs.high.y, cs.high.z);
	//}
	for (int s = 0; s < flexNBin; s++) {
		
		/*
		// test if the bin values are 0
		if (cornerSumHistogram[0][s] > 100000) { printf("cornerSumHistogram[0][%2d] = %f\n", s, cornerSumHistogram[0][s]); }
		if (cornerSumHistogram[1][s] > 100000) { printf("cornerSumHistogram[1][%2d] = %f\n", s, cornerSumHistogram[1][s]); }
		if (cornerSumHistogram[2][s] > 100000) { printf("cornerSumHistogram[2][%2d] = %f\n", s, cornerSumHistogram[2][s]); }
		if (cornerSumHistogram[3][s] > 100000) { printf("cornerSumHistogram[3][%2d] = %f\n", s, cornerSumHistogram[3][s]); }
		if (cornerSumHistogram[4][s] > 100000) { printf("cornerSumHistogram[4][%2d] = %f\n", s, cornerSumHistogram[4][s]); }
		if (cornerSumHistogram[5][s] > 100000) { printf("cornerSumHistogram[5][%2d] = %f\n", s, cornerSumHistogram[5][s]); }
		if (cornerSumHistogram[6][s] > 100000) { printf("cornerSumHistogram[6][%2d] = %f\n", s, cornerSumHistogram[6][s]); }
		if (cornerSumHistogram[7][s] > 100000) { printf("cornerSumHistogram[7][%2d] = %f\n", s, cornerSumHistogram[7][s]); }
		*/

		blockHistogram[s] = cornerSumHistogram[0][s] + cornerSumHistogram[5][s] + cornerSumHistogram[2][s] + cornerSumHistogram[7][s] - cornerSumHistogram[1][s] - cornerSumHistogram[4][s] - cornerSumHistogram[3][s] - cornerSumHistogram[6][s]; 
		if (blockHistogram[s] < 0) { blockHistogram[s] = 0; }	// NOTE: here we just make sure it has no negative value
		if (n == 5 && blockHistogram[s] > 0.000001) {
			printf("i = %2d, freq = %1.6f\n", s, blockHistogram[s]);
		}
	}

	// normalize blockHistogram
	Span currentSpan = flexBlock[n];
	//int weight = (currentSpan.high.x - currentSpan.low.x + 1) * (currentSpan.high.y - currentSpan.low.y + 1) * (currentSpan.high.z - currentSpan.low.z + 1);
	//if (n == 5) {
	//	printf("n = %2d, spanLow(%2d, %2d, %2d), spanHigh(%2d, %2d, %2d), weight = %d\n", n, currentSpan.low.x, currentSpan.low.y, currentSpan.low.z, currentSpan.high.x, currentSpan.high.y, currentSpan.high.z, weight);
	//}
	float totalBlockHistogram = 0;
	for (int s = 0; s < flexNBin; s++) {
		totalBlockHistogram += blockHistogram[s];
	}
	if (totalBlockHistogram <= 0) {
		printf("block %d: totalBlockHistogram = %f\n", n, totalBlockHistogram);
	}
	else {
		for (int s = 0; s < flexNBin; s++) {
			blockHistogram[s] = blockHistogram[s] / totalBlockHistogram;
			if (blockHistogram[s] < 0) { blockHistogram[s] = 0; }
			if (blockHistogram[s] > 1) { blockHistogram[s] = 1; }
			//if (n >= 24 && n < 27) {
			//	printf("%1.6f\n", blockHistogram[s]);
			//}
		}
	}


	//compute entropy
	float entropy = 0;
	for( int i = 0; i < flexNBin; i++ )
	{
		float probability = blockHistogram[i];
		//printf("%2d %d %f\t", n, i, probability);
		entropy += ( probability * ( probability <= 0 ? 0 : ( log( probability ) / log(2.0) ) ) );
	}

	// Change sign
	entropy = -entropy;

	// Normalize, if required
	entropy /= ( log( (float)flexNBin ) / log( 2.0f ) );

	flexBlockData[n].z = entropy;
	//printf("blockId = %2d, entropy = %f\n", n, entropy);
}

// d_queryBlockNew is to get subSpan of each corner
__global__ void
d_queryBlockNew(int blockNumber) {
	//printf("blockIdx: %d, threadIdx: %d\n", blockIdx.x, threadIdx.x);
	//int4 SpanLow = make_int4(15, 20, 25, 0); 
	//int4 SpanHigh = make_int4(21, 28, 33, 0);
	int4 SpanLow = flexBlock[blockNumber].low;
	int4 SpanHigh = flexBlock[blockNumber].high;
	int cornerId = blockIdx.x;
	//int spanId = threadIdx.x;
	switch (cornerId) {
	// 2: (0, 1, 0)		3: (1, 1, 0)
	// 0: (0, 0, 0)		1: (1, 0, 0)
	// 6: (0, 1, 1)		7: (1, 1, 1)
	// 4: (0, 0, 1)		5: (1, 0, 1)
	case 0:
		corner[0].x = (SpanLow.x < SpanHigh.x) ? SpanLow.x : SpanHigh.x;
		corner[0].y = (SpanLow.y < SpanHigh.y) ? SpanLow.y : SpanHigh.y;
		corner[0].z = (SpanLow.z < SpanHigh.z) ? SpanLow.z : SpanHigh.z;
		corner[0].w = 0;
		break;
	case 1:
		corner[1].x = (SpanLow.x > SpanHigh.x) ? SpanLow.x : SpanHigh.x;
		corner[1].y = (SpanLow.y < SpanHigh.y) ? SpanLow.y : SpanHigh.y;
		corner[1].z = (SpanLow.z < SpanHigh.z) ? SpanLow.z : SpanHigh.z;
		corner[1].w = 0;
		break;
	case 2:
		corner[2].x = (SpanLow.x < SpanHigh.x) ? SpanLow.x : SpanHigh.x;
		corner[2].y = (SpanLow.y > SpanHigh.y) ? SpanLow.y : SpanHigh.y;
		corner[2].z = (SpanLow.z < SpanHigh.z) ? SpanLow.z : SpanHigh.z;
		corner[2].w = 0;
		break;
	case 3:
		corner[3].x = (SpanLow.x > SpanHigh.x) ? SpanLow.x : SpanHigh.x;
		corner[3].y = (SpanLow.y > SpanHigh.y) ? SpanLow.y : SpanHigh.y;
		corner[3].z = (SpanLow.z < SpanHigh.z) ? SpanLow.z : SpanHigh.z;
		corner[3].w = 0;
		break;
	case 4:
		corner[4].x = (SpanLow.x < SpanHigh.x) ? SpanLow.x : SpanHigh.x;
		corner[4].y = (SpanLow.y < SpanHigh.y) ? SpanLow.y : SpanHigh.y;
		corner[4].z = (SpanLow.z > SpanHigh.z) ? SpanLow.z : SpanHigh.z;
		corner[4].w = 0;
		break;
	case 5:
		corner[5].x = (SpanLow.x > SpanHigh.x) ? SpanLow.x : SpanHigh.x;
		corner[5].y = (SpanLow.y < SpanHigh.y) ? SpanLow.y : SpanHigh.y;
		corner[5].z = (SpanLow.z > SpanHigh.z) ? SpanLow.z : SpanHigh.z;
		corner[5].w = 0;
		break;
	case 6:
		corner[6].x = (SpanLow.x < SpanHigh.x) ? SpanLow.x : SpanHigh.x;
		corner[6].y = (SpanLow.y > SpanHigh.y) ? SpanLow.y : SpanHigh.y;
		corner[6].z = (SpanLow.z > SpanHigh.z) ? SpanLow.z : SpanHigh.z;
		corner[6].w = 0;
		break;
	case 7:
		corner[7].x = (SpanLow.x > SpanHigh.x) ? SpanLow.x : SpanHigh.x;
		corner[7].y = (SpanLow.y > SpanHigh.y) ? SpanLow.y : SpanHigh.y;
		corner[7].z = (SpanLow.z > SpanHigh.z) ? SpanLow.z : SpanHigh.z;
		corner[7].w = 0;
		break;
	}
	//if (threadIdx.x == 0) {
	//	printf("corner[%d] %d %d %d %d\n", cornerId, corner[cornerId].x, corner[cornerId].y, corner[cornerId].z, corner[cornerId].w);
	//}
	int2 subSpanX[6];	// NOTE: here we hard code array size 6, because bin number is 64 = 2^6,
					// there are at most 6 sub spans in each dimension
					// subSpanX[i].x is the lower span, subSpanX[i].y is the upper span
	int2 subSpanY[6];
	int2 subSpanZ[6];
	//Span subSpan[6*6*6]; // at most 6*6*6 sub spans
	int nx = 0;			// count how many spans generated in each dimension
	int ny = 0;
	int nz = 0;
	int n = 0;	// total number of subSpans
	int x = corner[cornerId].x;
	int y = corner[cornerId].y;
	int z = corner[cornerId].z;
	//printf("corner[%d] %d %d %d\n", t, corners[t].x, corners[t].y, corners[t].z);
		
	for (int i = 0; i < 6; i++) {
		if ((x & (~(1 << i))) != x) {
			subSpanX[nx].y = x;
			x &= ~(1 << i);
			subSpanX[nx].x = x + 1;
			//printf("subSpanX[%d]: low %d, high %d\n", nx, subSpanX[nx].x, subSpanX[nx].y);
			nx++;
		}
		if (x == 0) { break; }
	}
	for (int i = 0; i < 6; i++) {
		if ((y & (~(1 << i))) != y) {
			subSpanY[ny].y = y;
			y &= ~(1 << i);
			subSpanY[ny].x = y + 1;
			//printf("subSpanY[%d]: low %d, high %d\n", ny, subSpanY[ny].x, subSpanY[ny].y);
			ny++;
		}
		if (y == 0) { break; }
	}
	for (int i = 0; i < 6; i++) {
		if ((z & (~(1 << i))) != z) {
			subSpanZ[nz].y = z;
			z &= ~(1 << i);
			subSpanZ[nz].x = z + 1;
			//printf("subSpanZ[%d]: low %d, high %d\n", nz, subSpanZ[nz].x, subSpanZ[nz].y);
			nz++;
		}
		if (z == 0) { break; }
	}
	n = nx * ny * nz;
	nSubSpan[cornerId] = n;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				subSpan[cornerId][i*ny*nz + j*nz + k].low.x = subSpanX[i].x;	
				subSpan[cornerId][i*ny*nz + j*nz + k].low.y = subSpanY[j].x;
				subSpan[cornerId][i*ny*nz + j*nz + k].low.z = subSpanZ[k].x;
				subSpan[cornerId][i*ny*nz + j*nz + k].high.x = subSpanX[i].y;
				subSpan[cornerId][i*ny*nz + j*nz + k].high.y = subSpanY[j].y;
				subSpan[cornerId][i*ny*nz + j*nz + k].high.z = subSpanZ[k].y;
			}
		}
	}
	
	if (threadIdx.x == 0) {
		/*
		for (int i = 0; i < n; i++) {
			printf("subSpan[%2d]: low(%2d, %2d, %2d)  high(%2d, %2d, %2d)\n", i,
					subSpan[cornerId][i].low.x, subSpan[cornerId][i].low.y, subSpan[cornerId][i].low.z,
					subSpan[cornerId][i].high.x, subSpan[cornerId][i].high.y, subSpan[cornerId][i].high.z);
		}
		*/
	//	printf("nSubSpan[%d] (%2d, %2d, %2d) n = %d\n", cornerId, corner[cornerId].x, corner[cornerId].y, corner[cornerId].z, nSubSpan[cornerId]);
	}


	//float subSpanHistogram[flexNBin];	// this array stores the sum of all subSpan histograms of each corner
	//for (int countSSH = 0; countSSH < flexNBin; countSSH++) {
	//	subSpanHistogram[countSSH] = 0;
	//}
	// look up histogram of n subSpans
	//printf("spanLow %d %d %d, spanHigh %d %d %d, ", SpanLow.x, SpanLow.y, SpanLow.z, SpanHigh.x, SpanHigh.y, SpanHigh.z);
	//printf("current t %d, n %d, nx %d, ny %d, nz %d\n", t, n, nx, ny, nz);
}

// d_querySpanNew is to get histogram of each span
__global__ void
d_querySpanNew() {

	int cId = blockIdx.x;	// corner Id
	int sId = threadIdx.x;	// span Id
	if (sId >= nSubSpan[cId])
		return;

	Span currentSpan = subSpan[cId][sId];
	int4 currentFractal;	// current fractal encoded histogram
	int indexError = 0;		// index of NE in errorTexture
			
	// NOTE: Here we need to use a threshold to distinguish whether this span is fractal encoded or simple histogram
	//if (d_spanSize(currentSpan) >= 8) {
	if (1) 
	{	
		// look up fractal codebook
		int found = 0;
		for (int iz = 0; iz < nDimension/2; iz++) {
			for (int iy = 0; iy < nDimension; iy++) {
				for (int ix = 0; ix < nDimension; ix++) {	// NOTE: the fuel dataset size is 64x64x64, but half of them are encoded by fracal and half by simple histogram, so that each texture has 64x64x32 size
					int4 spanLow = tex3D(codebookSpanLowTex, ix, iy, iz);
					int4 spanHigh = tex3D(codebookSpanHighTex, ix, iy, iz);		
					if ((spanLow.x == currentSpan.low.x)
							&& (spanLow.y == currentSpan.low.y)
							&& (spanLow.z == currentSpan.low.z)
							&& (spanHigh.x == currentSpan.high.x)
							&& (spanHigh.y == currentSpan.high.y)
							&& (spanHigh.z == currentSpan.high.z)) {
						currentFractal = tex3D(flexibleCodebookTex, ix, iy, iz);
						// TODO: add NE
						//indexError = ix * iy * iz;
						// NOTE: how to index error list?
						indexError = ix + iy * nDimension + iz * nDimension * nDimension;
						//printf("cId %d, sId %2d, currentFractal: ID %3d, shift %2d, flip %d, NE %d\n", cId, sId, currentFractal.x, currentFractal.y, currentFractal.z, currentFractal.w);
						found = 1;
						break;
					}
				}
			}
		}// end look up fractal codebook
		if (found == 0) {
			printf("didn't find fractal: cId %d, sId %2d, spanLow(%2d, %2d, %2d), spanHigh(%2d, %2d, %2d)\n", cId, sId, currentSpan.low.x, currentSpan.low.y, currentSpan.low.z, currentSpan.high.x, currentSpan.high.y, currentSpan.high.z);
		}
		
		// fractal decoding
		// Flip (if needed)->Shift->merge errors->normalize the frequencies
		float originalTemplate[flexNBin];
		float currentTemplate[flexNBin];
		int templateId = currentFractal.x;
		// TODO: find out why next line makes program crush
		if (templateId < 0 || templateId > nTemplate) { printf("Error! templateID: %d\n", templateId); }
		int shift = currentFractal.y;
		int flipFlag = currentFractal.z;
		int NE = currentFractal.w;
		// TODO: find out why we cannot print NE
		if (NE < 0 || NE > flexNBin) {
			printf("Error NE!: %d\n", NE);
		}
		//float sumOT = 0;
		for (int i = 0; i < flexNBin; i++) {
		// TODO: find out why next line makes program crush
			originalTemplate[i] = tex2DLayered(flexibleTemplatesTex, (float) i, (float) templateId, 1);
			if (originalTemplate[i] < 0 || originalTemplate[i] > 1) {
				printf("Error! originalTemplate[%d]: %f\n", i, originalTemplate[i]);
			}
			//sumOT += originalTemplate[i];
		}
		//printf("cId %d, sId %2d, sumOT %f\n", cId, sId, sumOT);
		flexibleFractalDecoding( originalTemplate, currentTemplate, flexNBin, flipFlag, shift);
		for (int i = 0; i < flexNBin; i++) {
			if (currentTemplate[i] < 0 || currentTemplate[i] > 1) {
				printf("Error! currentTemplate[%d], %f\n", i, currentTemplate[i]);
			}
		}
		// add NE
		for (int i = 0; i < NE; i++) {
			// TODO: find out why we could not look up texture and use if printf command
			float2 error = tex2DLayered(flexibleErrorsbookTex, (float) i, (float)(indexError % 2048), (int) (indexError / 2048) );
			int errorIndex = (int)error.x;
			if (errorIndex < 0 || errorIndex > flexNBin) {
				printf("Error Index! %d\n", errorIndex);
			}
			float errorValue = error.y;
			if (errorValue < -1 || errorValue > 1) {
				printf("Error Value! %f\n", errorValue);
			}
			// TODO: find out why next line makes program crush
			currentTemplate[errorIndex] += errorValue;	// TODO: in this step, we get wrong value
			if (currentTemplate[errorIndex] < 0 ) { currentTemplate[errorIndex] = 0; }
		}
		// normalize currenTemplate
		float tempTotal = 0;
		for (int i = 0; i < flexNBin; i++) {
			tempTotal += currentTemplate[i];
		}
		for (int i = 0; i < flexNBin; i++) {
			currentTemplate[i] = currentTemplate[i] / tempTotal;
			if (currentTemplate[i] < 0 || currentTemplate[i] > 1) { 
				printf("Error! currentTemplate[%d]: %f\n", i, currentTemplate[i]); 
			}
		}

		/*
		// test currentTemplate
		float sumCT = 0;
		for (int i = 0; i < flexNBin; i++) {
			sumCT += currentTemplate[i];
		}
		if (sumCT > 1.000001 || sumCT < 0.999999) {
			printf("cId %d, sId %2d, sumCT: %f\n", cId, sId, sumCT);
		}*/
		

		// copy currentTemplate into cornerHistogram
		for (int countCopy = 0; countCopy < flexNBin; countCopy++) {
			cornerHistogram[cId][sId][countCopy] = currentTemplate[countCopy];
		}

		
		// test currentTemplate
		float sumCT = 0;
		for (int i = 0; i < flexNBin; i++) {
			sumCT += cornerHistogram[cId][sId][i];
		}
		if (sumCT > 1.000001 || sumCT < 0.999999) {
			printf("cId %d, sId %2d, sumCT: %f\n", cId, sId, sumCT);
		}

		// NOTE: we could not directly copy currentTemplate to cornerHistogram, beacuse they need synchronize
		// TODO: use shared memory
	}
	else {
		//look up simple codebook
		// look up simple histogram
		//printf("simple! cId %d, sId %2d, low(%2d, %2d, %2d), high(%2d, %2d, %2d)\n", cId, sId, currentSpan.low.x, currentSpan.low.y, currentSpan.low.z, currentSpan.high.x, currentSpan.high.y, currentSpan.high.z);
		
		// NOTE: when we are doing bitwise operation, we assume range (1 ~ 64), but the simple span data range is (0 ~ 63)
		// NOTE: but the fractal span range is still (1 ~ 64)
		currentSpan.low.x -= 1;
		currentSpan.low.y -= 1;
		currentSpan.low.z -= 1;
		currentSpan.high.x -= 1;
		currentSpan.high.y -= 1;
		currentSpan.high.z -= 1;
		
		int currentSimpleCount = 0;
		int indexSimpleHistogram = 0;

		float currentSimpleHistogram[flexNBin];
		// initialize
		for (int i = 0; i < flexNBin; i++) {
			currentSimpleHistogram[i] = 0;
		}
		
		int found = 0;
		
		for (int iz = 0; iz < nDimension/2; iz++) {
			for (int iy = 0; iy < nDimension; iy++) {
				for (int ix = 0; ix < nDimension; ix++) {
					int4 spanLow = tex3D(simpleSpanLowTex, ix, iy, iz);
					int4 spanHigh = tex3D(simpleSpanHighTex, ix, iy, iz);
					//if (ix < 5 && iy < 5 && iz < 5) {
					//	printf("%d %d %d, spanLow(%d, %d, %d) spanHigh(%d, %d, %d)\n", ix, iy, iz, spanLow.x, spanLow.y, spanLow.z, spanHigh.x, spanHigh.y, spanHigh.z);
					//}
					if (spanLow.x == currentSpan.low.x
							&& spanLow.y == currentSpan.low.y
							&& spanLow.z == currentSpan.low.z
							&& spanHigh.x == currentSpan.high.x
							&& spanHigh.y == currentSpan.high.y
							&& spanHigh.z == currentSpan.high.z) {
						currentSimpleCount = tex3D(simpleCountTex, ix, iy, iz);
						//indexSimpleHistogram = ix * iy * iz;
						// NOTE: how to index simple histogram?
						indexSimpleHistogram = ix + iy * nDimension + iz * nDimension * nDimension;
						//printf("cId %d, sId %2d, Count %2d, index %d\n", cId, sId, currentSimpleCount, indexSimpleHistogram);
						found = 1;
						break;
					}
				}
			}
		}// end looking up simpleHistogram
		if (found == 0) {
			printf("didn't find simple! cId %d, sId %2d, simpleLow(%2d, %2d, %2d), simpleHigh(%2d, %2d, %2d)\n", cId, sId, currentSpan.low.x, currentSpan.low.y, currentSpan.low.z, currentSpan.high.x, currentSpan.high.y, currentSpan.high.z);
		}
		for (int ih = 0; ih < currentSimpleCount; ih++) {
			float2 current = tex2DLayered(simpleHistogramTex, ih, indexSimpleHistogram % 2048, indexSimpleHistogram / 2048);
			currentSimpleHistogram[(int)current.x] = current.y;
			//printf("cId %d, sId %d, bin %d, freq %f\n", cId, sId, current.x, current.y);
		}
		for (int countCopy = 0; countCopy < flexNBin; countCopy++) {
			cornerHistogram[cId][sId][countCopy] = currentSimpleHistogram[countCopy];
		}
		
		
		// test currentTemplate
		float sumCT = 0;
		for (int i = 0; i < flexNBin; i++) {
			sumCT += cornerHistogram[cId][sId][i];
		}
		if (sumCT > 1.000001 || sumCT < 0.999999) {
			printf("cId %d, sId %2d, sumCT: %f\n", cId, sId, sumCT);
		}
		
	}// end look up simple codebook
}

__global__ void
d_clearCornerHistogram() {
	int cId = blockIdx.x;
	int tId = threadIdx.x;
//	cornerHistogram[cId][tId] = 0;
}

void dataProcessing() {
	//d_queryBlock<<<1, 1>>>(make_int4(15, 20, 25, 0), make_int4(15, 20, 25, 0));
	//d_queryBlockNew<<<8, 1>>>(make_int4(15, 20, 25, 0), make_int4(21, 28, 33, 0));
	//d_queryBlockNew<<<8, 1>>>(3);
	//d_querySpanNew<<<8, 216>>>();
	d_divideBlock<<<1, 1>>>(30);
	int h_nFlexBlock;	// how many blocks in this volume
	cudaMemcpyFromSymbol(&h_nFlexBlock, nFlexBlock, sizeof(int));
	printf("h_nFlexBlock: %d\n", h_nFlexBlock);
	
	for (int i = 0; i < h_nFlexBlock; i++) {
		//d_clearCornerHistogram<<<8, flexNBin>>>();
		d_queryBlockNew<<<8, 1>>>(i);
		d_querySpanNew<<<8, 216>>>();
		d_sumSanHistogram<<<8, 1>>>(i);
		//d_normalizeCorner<<<8, 216>>>();
		d_computeBlock<<<1, 1>>>(i);
		//TODO: store the result into an array, and pass that array into a texture
	}
	//d_testSimpleSum<<<1, 1>>>();
	cudaDeviceSynchronize();
	//d_testEntropy<<<1, 1>>>();
	//bindToTex();
	//d_testTexRef<<<1, 1>>>();
	//bindToTexObj();
	//d_testTexObj<<<1, 1>>>(flexBlockTexObj);

	//bindToSurfObj();
	//d_testSurfObj<<<1, 1>>>(flexBlockSurfObj);

	//d_testMallocMem<<<1, 1>>>();
	//d_testTd<<<3, 3>>>();
	//d_testSharedMemory<<<1, 4>>>();
	//d_testPrintShared<<<1, 1>>>();
	//d_showEntropy<<<1, h_nFlexBlock>>>();
	//d_testSimpleSpan<<<64, 64>>>();
}


void 
initCuda( int4 *h_codebookSpanLow, int4 *h_codebookSpanHigh,
		  int4 *h_flexibleCodebook, 
		  float2 *h_flexibleErrorsbook, 
		  int4 *h_simpleLow, int4 *h_simpleHigh, int *h_simpleCount,
		  float2 *h_simpleHistogram, float *h_flexibleTemplates )
{
	// create a volume data in which each voxel has a index number
	// which can be used to lookup the histogram location
	//int volume[nBlocks];
	//for (int i = 0; i < nBlocks; i++)
	//	volume[i] = i;
	
    ////////////////////////////////////////////
    // flexible block size
    ////////////////////////////////////////////

    // create 3D array for codebookSpanLow
	cudaChannelFormatDesc channelDescCodebookSpanLow = cudaCreateChannelDesc<SpanType>();
    checkCudaErrors(cudaMalloc3DArray(&d_codebookSpanLowArray, &channelDescCodebookSpanLow, flexibleVolumeSize));

    // create 3D array for codebookSpanHigh
    cudaChannelFormatDesc channelDescCodebookSpanHigh = cudaCreateChannelDesc<SpanType>();
    checkCudaErrors(cudaMalloc3DArray(&d_codebookSpanHighArray, &channelDescCodebookSpanHigh, flexibleVolumeSize));

    // create 3D array for flexibleCodebook
    cudaChannelFormatDesc channelDescFlexibleCodebook = cudaCreateChannelDesc<FlexibleCodebookType>();
    checkCudaErrors(cudaMalloc3DArray(&d_flexibleCodebookArray, &channelDescFlexibleCodebook, flexibleVolumeSize));

    // create 3D array for flexibleErrorsbook
    cudaChannelFormatDesc channelDescFlexibleErrorsbook = cudaCreateChannelDesc<FlexibleErrorsbookType>();
    checkCudaErrors(cudaMalloc3DArray(&d_flexibleErrorsbookArray, &channelDescFlexibleErrorsbook, flexibleHistogramSize, cudaArrayLayered));

    // create 3D array for simpleSpanLow
    cudaChannelFormatDesc channelDescSimpleSpanLow = cudaCreateChannelDesc<SpanType>();
    checkCudaErrors(cudaMalloc3DArray(&d_simpleSpanLowArray, &channelDescSimpleSpanLow, flexibleVolumeSize));

    // create 3D array for simpleSpanHigh
    cudaChannelFormatDesc channelDescSimpleSpanHigh = cudaCreateChannelDesc<SpanType>();
    checkCudaErrors(cudaMalloc3DArray(&d_simpleSpanHighArray, &channelDescSimpleSpanHigh, flexibleVolumeSize));

    // create 3D array for simpleCount
    cudaChannelFormatDesc channelDescSimpleCount = cudaCreateChannelDesc<SimpleCountType>();
    checkCudaErrors(cudaMalloc3DArray(&d_simpleCountArray, &channelDescSimpleCount, flexibleVolumeSize));

    // create 3D array for simpleHistogram
    cudaChannelFormatDesc channelDescSimpleHistogram = cudaCreateChannelDesc<SimpleHistogramType>();
    checkCudaErrors(cudaMalloc3DArray(&d_simpleHistogramArray, &channelDescSimpleHistogram, flexibleHistogramSize, cudaArrayLayered));

    // create 3D array for flexibleTemplates
	cudaChannelFormatDesc channelDescFlexibleTemplates = cudaCreateChannelDesc<TemplatesType>();
    checkCudaErrors(cudaMalloc3DArray(&d_flexibleTemplatesArray, &channelDescFlexibleTemplates, flexibleTemplatesSize, cudaArrayLayered));

	////////////////////////////////////////////////
	// flexible block size
	///////////////////////////////////////////////

	// copy data to 3D array of codebookSpanLow
	cudaMemcpy3DParms copyParamsCodebookSpanLow = {0};
	copyParamsCodebookSpanLow.srcPtr	= make_cudaPitchedPtr(h_codebookSpanLow, flexibleVolumeSize.width*sizeof(SpanType), flexibleVolumeSize.width, flexibleVolumeSize.height);
	copyParamsCodebookSpanLow.dstArray	= d_codebookSpanLowArray;
	copyParamsCodebookSpanLow.extent	= flexibleVolumeSize;
	copyParamsCodebookSpanLow.kind		= cudaMemcpyHostToDevice;
	checkCudaErrors(cudaMemcpy3D(&copyParamsCodebookSpanLow));

	// copy data to 3D array of codebookSpanHigh
	cudaMemcpy3DParms copyParamsCodebookSpanHigh = {0};
	copyParamsCodebookSpanHigh.srcPtr	= make_cudaPitchedPtr(h_codebookSpanHigh, flexibleVolumeSize.width*sizeof(SpanType), flexibleVolumeSize.width, flexibleVolumeSize.height);
	copyParamsCodebookSpanHigh.dstArray	= d_codebookSpanHighArray;
	copyParamsCodebookSpanHigh.extent	= flexibleVolumeSize;
	copyParamsCodebookSpanHigh.kind		= cudaMemcpyHostToDevice;
	checkCudaErrors(cudaMemcpy3D(&copyParamsCodebookSpanHigh));

	// copy data to 3D array of flexibleCodebook
	cudaMemcpy3DParms copyParamsFlexibleCodebook = {0};
	copyParamsFlexibleCodebook.srcPtr	= make_cudaPitchedPtr(h_flexibleCodebook, flexibleVolumeSize.width*sizeof(FlexibleCodebookType), flexibleVolumeSize.width, flexibleVolumeSize.height);
	copyParamsFlexibleCodebook.dstArray	= d_flexibleCodebookArray;
	copyParamsFlexibleCodebook.extent	= flexibleVolumeSize;
	copyParamsFlexibleCodebook.kind		= cudaMemcpyHostToDevice;
	checkCudaErrors(cudaMemcpy3D(&copyParamsFlexibleCodebook));

	// copy data to 3D array of flexibleErrorsbook
	cudaMemcpy3DParms copyParamsFlexibleErrorsbook = {0};
	copyParamsFlexibleErrorsbook.srcPtr 	= make_cudaPitchedPtr(h_flexibleErrorsbook, flexibleHistogramSize.width*sizeof(FlexibleErrorsbookType), flexibleHistogramSize.width, flexibleHistogramSize.height);
	copyParamsFlexibleErrorsbook.dstArray	= d_flexibleErrorsbookArray;
	copyParamsFlexibleErrorsbook.extent		= flexibleHistogramSize;
	copyParamsFlexibleErrorsbook.kind		= cudaMemcpyHostToDevice;

	checkCudaErrors(cudaMemcpy3D(&copyParamsFlexibleErrorsbook));

	// copy data to 3D array of simpleLow
	cudaMemcpy3DParms copyParamsSimpleLow = {0};
	copyParamsSimpleLow.srcPtr		= make_cudaPitchedPtr(h_simpleLow, flexibleVolumeSize.width*sizeof(SpanType), flexibleVolumeSize.width, flexibleVolumeSize.height);
	copyParamsSimpleLow.dstArray	= d_simpleSpanLowArray;
	copyParamsSimpleLow.extent		= flexibleVolumeSize;
	copyParamsSimpleLow.kind		= cudaMemcpyHostToDevice;
	checkCudaErrors(cudaMemcpy3D(&copyParamsSimpleLow));

	// copy data to 3D array of simpleHigh
	cudaMemcpy3DParms copyParamsSimpleHigh = {0};
	copyParamsSimpleHigh.srcPtr		= make_cudaPitchedPtr(h_simpleHigh, flexibleVolumeSize.width*sizeof(SpanType), flexibleVolumeSize.width, flexibleVolumeSize.height);
	copyParamsSimpleHigh.dstArray	= d_simpleSpanHighArray;
	copyParamsSimpleHigh.extent		= flexibleVolumeSize;
	copyParamsSimpleHigh.kind		= cudaMemcpyHostToDevice;
	checkCudaErrors(cudaMemcpy3D(&copyParamsSimpleHigh));

	// copy data to 3D array of simpleCount
	cudaMemcpy3DParms copyParamsSimpleCount = {0};
	copyParamsSimpleCount.srcPtr	= make_cudaPitchedPtr(h_simpleCount, flexibleVolumeSize.width*sizeof(SimpleCountType), flexibleVolumeSize.width, flexibleVolumeSize.height);
	copyParamsSimpleCount.dstArray	= d_simpleCountArray;
	copyParamsSimpleCount.extent	= flexibleVolumeSize;
	copyParamsSimpleCount.kind		= cudaMemcpyHostToDevice;
	checkCudaErrors(cudaMemcpy3D(&copyParamsSimpleCount));

	// copy data to 3D array of simpleHistogram
	cudaMemcpy3DParms copyParamsSimpleHistogram = {0};
	copyParamsSimpleCount.srcPtr	= make_cudaPitchedPtr(h_simpleHistogram, flexibleHistogramSize.width*sizeof(SimpleHistogramType), flexibleHistogramSize.width, flexibleHistogramSize.height);
	copyParamsSimpleCount.dstArray	= d_simpleHistogramArray;
	copyParamsSimpleCount.extent	= flexibleHistogramSize;
	copyParamsSimpleCount.kind		= cudaMemcpyHostToDevice;
	checkCudaErrors(cudaMemcpy3D(&copyParamsSimpleCount));

	// copy data to 3D array of flexibleTemplates
	cudaMemcpy3DParms copyParamsFlexibleTemplates = {0};
	copyParamsFlexibleTemplates.srcPtr		= make_cudaPitchedPtr(h_flexibleTemplates, flexibleTemplatesSize.width*sizeof(TemplatesType), flexibleTemplatesSize.width, flexibleTemplatesSize.height);
	copyParamsFlexibleTemplates.dstArray	= d_flexibleTemplatesArray;
	copyParamsFlexibleTemplates.extent		= flexibleTemplatesSize;
	copyParamsFlexibleTemplates.kind		= cudaMemcpyHostToDevice;
	checkCudaErrors(cudaMemcpy3D(&copyParamsFlexibleTemplates));

	///////////////////////////////////////////////
	// flexible block size
	///////////////////////////////////////////////

	// set texture parameters for codebookSpanLowTex
	codebookSpanLowTex.normalized = false;
	codebookSpanLowTex.filterMode = cudaFilterModePoint;
	codebookSpanLowTex.addressMode[0] = cudaAddressModeClamp;
	codebookSpanLowTex.addressMode[1] = cudaAddressModeClamp;

	// set texture parameters for codebookSpanHighTex
	codebookSpanHighTex.normalized = false;
	codebookSpanHighTex.filterMode = cudaFilterModePoint;
	codebookSpanHighTex.addressMode[0] = cudaAddressModeClamp;
	codebookSpanHighTex.addressMode[1] = cudaAddressModeClamp;

	// set texture parameters for flexibleCodebookTex
	flexibleCodebookTex.normalized = false;
	flexibleCodebookTex.filterMode = cudaFilterModePoint;
	flexibleCodebookTex.addressMode[0] = cudaAddressModeClamp;
	flexibleCodebookTex.addressMode[1] = cudaAddressModeClamp;

	// set texture parameters for flexibleErrorsbookTex
	flexibleErrorsbookTex.normalized = false;
	flexibleErrorsbookTex.filterMode = cudaFilterModePoint;
	flexibleErrorsbookTex.addressMode[0] = cudaAddressModeClamp;
	flexibleErrorsbookTex.addressMode[1] = cudaAddressModeClamp;

	// set texture parameters for simpleSpanLowTex
	simpleSpanLowTex.normalized = false;
	simpleSpanLowTex.filterMode = cudaFilterModePoint;
	simpleSpanLowTex.addressMode[0] = cudaAddressModeClamp;
	simpleSpanLowTex.addressMode[1] = cudaAddressModeClamp;

	// set texture parameters for simpleSpanHighTex
	simpleSpanHighTex.normalized = false;
	simpleSpanHighTex.filterMode = cudaFilterModePoint;
	simpleSpanHighTex.addressMode[0] = cudaAddressModeClamp;
	simpleSpanHighTex.addressMode[1] = cudaAddressModeClamp;

	// set texture parameters for simpleCountTex
	simpleCountTex.normalized = false;
	simpleCountTex.filterMode = cudaFilterModePoint;
	simpleCountTex.addressMode[0] = cudaAddressModeClamp;
	simpleCountTex.addressMode[1] = cudaAddressModeClamp;

	// set texture parameters for simpleHistogramTex
	simpleHistogramTex.normalized = false;
	simpleHistogramTex.filterMode = cudaFilterModePoint;
	simpleHistogramTex.addressMode[0] = cudaAddressModeClamp;
	simpleHistogramTex.addressMode[1] = cudaAddressModeClamp;

	// set texture parameters for flexibleTemplatesTex
	flexibleTemplatesTex.normalized = false;
	flexibleTemplatesTex.filterMode = cudaFilterModePoint;
	flexibleTemplatesTex.addressMode[0] = cudaAddressModeClamp;
	flexibleTemplatesTex.addressMode[1] = cudaAddressModeClamp;

	/////////////////////////////////////////////////////////////
	// flexible block size
	/////////////////////////////////////////////////////////////

	// bind array to 3D texture of codebookSpanLow
	checkCudaErrors(cudaBindTextureToArray(codebookSpanLowTex, d_codebookSpanLowArray, channelDescCodebookSpanLow));

	// bind array to 3D texture of codebookSpanHigh
	checkCudaErrors(cudaBindTextureToArray(codebookSpanHighTex, d_codebookSpanHighArray, channelDescCodebookSpanHigh));

	// bind array to 3D texture of flexibleCodebook
	checkCudaErrors(cudaBindTextureToArray(flexibleCodebookTex, d_flexibleCodebookArray, channelDescFlexibleCodebook));

	// bind array to 2D layered texture of flexibleErrorsbook
	checkCudaErrors(cudaBindTextureToArray(flexibleErrorsbookTex, d_flexibleErrorsbookArray, channelDescFlexibleErrorsbook));

	// bind array to 3D texture of simpleSpanLow
	checkCudaErrors(cudaBindTextureToArray(simpleSpanLowTex, d_simpleSpanLowArray, channelDescSimpleSpanLow));

	// bind array to 3D texture of simpleSpanHigh
	checkCudaErrors(cudaBindTextureToArray(simpleSpanHighTex, d_simpleSpanHighArray, channelDescSimpleSpanHigh));

	// bind array to 3D texture of simpleCount
	checkCudaErrors(cudaBindTextureToArray(simpleCountTex, d_simpleCountArray, channelDescSimpleCount));

	// bind array to 2D layered texture of simpleHistogram
	checkCudaErrors(cudaBindTextureToArray(simpleHistogramTex, d_simpleHistogramArray, channelDescSimpleHistogram));

	// bind array to 2D layered texture of flexibleTemplates
	checkCudaErrors(cudaBindTextureToArray(flexibleTemplatesTex, d_flexibleTemplatesArray, channelDescFlexibleTemplates));

	//for (int i = 0; i < 622; i++) {
	//	for (int j = 0; j < flexNBin; j++) {
	//		if (h_templates[i * flexNBin + j] < 0 || h_templates[i * flexNBin + j] > 1) {
	//			printf("ERROR! template(%d, %d): %f\n", i, j, h_templates[i * flexNBin + j]);
	//		}
	//		if (i == 32) {
	//			//printf("template(32, %d): %f\n", j, h_templates[i * nBins + j]);
	//		}
	//	}
	//}

}// end function

void freeCudaBuffers()
{
	// flexible block size
	checkCudaErrors(cudaFreeArray(d_codebookSpanLowArray));
	checkCudaErrors(cudaFreeArray(d_codebookSpanHighArray));
	checkCudaErrors(cudaFreeArray(d_flexibleCodebookArray));
	checkCudaErrors(cudaFreeArray(d_flexibleErrorsbookArray));
	checkCudaErrors(cudaFreeArray(d_simpleSpanLowArray));
	checkCudaErrors(cudaFreeArray(d_simpleSpanHighArray));
	checkCudaErrors(cudaFreeArray(d_simpleCountArray));
	checkCudaErrors(cudaFreeArray(d_simpleHistogramArray));
	checkCudaErrors(cudaFreeArray(d_flexibleTemplatesArray));

	checkCudaErrors(cudaFreeArray(d_flexTexArray));
}
