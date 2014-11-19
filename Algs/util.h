#ifndef SRMP_GHAODJFNADJGHAJSD
#define SRMP_GHAODJFNADJGHAJSD


// sorts arr[left..right]
template <class Item> inline void quickSort(Item* arr, int left, int right)
{
	while ( left < right )
	{
		int i = left, j = right;
		Item tmp;
		Item pivot = arr[(left + right) / 2];

		/* partition */
		while (i <= j)
		{
			while (arr[i] < pivot) i++;
			while (arr[j] > pivot) j--;
			if (i <= j)
			{
				tmp = arr[i];
				arr[i] = arr[j];
				arr[j] = tmp;
				i++;
				j--;
			}
		}
		/* recursion */
		if (2*j < left+right)
		{
			quickSort(arr, i, right);
			right = j;
		}
		else
		{
			quickSort(arr, left, j);
			left = i;
		}
	}
}

#endif
