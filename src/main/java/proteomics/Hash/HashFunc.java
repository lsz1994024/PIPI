package proteomics.Hash;

import java.nio.ByteBuffer;
import java.nio.charset.Charset;
import java.util.List;
import java.util.Map;

public class HashFunc {
    public static int hash32(String doc) {
        byte[] buffer = doc.getBytes(Charset.forName("utf-8"));
        ByteBuffer data = ByteBuffer.wrap(buffer);
        return hash32(data, 0, buffer.length, 0);
    }

    public static int hash32(ByteBuffer data, int offset, int length, int seed) {
        int m = 0x5bd1e995;
        int r = 24;

        int h = seed ^ length;

        int len_4 = length >> 2;

        for (int i = 0; i < len_4; i++) {
            int i_4 = i << 2;
            int k = data.get(offset + i_4 + 3);
            k = k << 8;
            k = k | (data.get(offset + i_4 + 2) & 0xff);
            k = k << 8;
            k = k | (data.get(offset + i_4 + 1) & 0xff);
            k = k << 8;
            k = k | (data.get(offset + i_4 + 0) & 0xff);
            k *= m;
            k ^= k >>> r;
            k *= m;
            h *= m;
            h ^= k;
        }

        // avoid calculating modulo
        int len_m = len_4 << 2;
        int left = length - len_m;

        if (left != 0) {
            if (left >= 3) {
                h ^= (int) data.get(offset + length - 3) << 16;
            }
            if (left >= 2) {
                h ^= (int) data.get(offset + length - 2) << 8;
            }
            if (left >= 1) {
                h ^= (int) data.get(offset + length - 1);
            }

            h *= m;
        }

        h ^= h >>> 13;
        h *= m;
        h ^= h >>> 15;

        return h;
    }

    public static long hash64(String doc) {
        byte[] buffer = doc.getBytes(Charset.forName("utf-8"));
        ByteBuffer data = ByteBuffer.wrap(buffer);
        return hash64(data, 0, buffer.length, 0);
    }

    public static long hash64(ByteBuffer key, int offset, int length, long seed) {
        long m64 = 0xc6a4a7935bd1e995L;
        int r64 = 47;

        long h64 = (seed & 0xffffffffL) ^ (m64 * length);

        int lenLongs = length >> 3;

        for (int i = 0; i < lenLongs; ++i) {
            int i_8 = i << 3;

            long k64 = ((long) key.get(offset + i_8 + 0) & 0xff)
                    + (((long) key.get(offset + i_8 + 1) & 0xff) << 8)
                    + (((long) key.get(offset + i_8 + 2) & 0xff) << 16)
                    + (((long) key.get(offset + i_8 + 3) & 0xff) << 24)
                    + (((long) key.get(offset + i_8 + 4) & 0xff) << 32)
                    + (((long) key.get(offset + i_8 + 5) & 0xff) << 40)
                    + (((long) key.get(offset + i_8 + 6) & 0xff) << 48)
                    + (((long) key.get(offset + i_8 + 7) & 0xff) << 56);

            k64 *= m64;
            k64 ^= k64 >>> r64;
            k64 *= m64;

            h64 ^= k64;
            h64 *= m64;
        }

        int rem = length & 0x7;

        switch (rem) {
            case 0:
                break;
            case 7:
                h64 ^= (long) key.get(offset + length - rem + 6) << 48;
            case 6:
                h64 ^= (long) key.get(offset + length - rem + 5) << 40;
            case 5:
                h64 ^= (long) key.get(offset + length - rem + 4) << 32;
            case 4:
                h64 ^= (long) key.get(offset + length - rem + 3) << 24;
            case 3:
                h64 ^= (long) key.get(offset + length - rem + 2) << 16;
            case 2:
                h64 ^= (long) key.get(offset + length - rem + 1) << 8;
            case 1:
                h64 ^= (long) key.get(offset + length - rem);
                h64 *= m64;
        }

        h64 ^= h64 >>> r64;
        h64 *= m64;
        h64 ^= h64 >>> r64;

        return h64;
    }

    public static int hammingDistance(int hash1, int hash2) {
        int i = hash1 ^ hash2;
        i = i - ((i >>> 1) & 0x55555555);
        i = (i & 0x33333333) + ((i >>> 2) & 0x33333333);
        i = (i + (i >>> 4)) & 0x0f0f0f0f;
        i = i + (i >>> 8);
        i = i + (i >>> 16);
        return i & 0x3f;
    }

    public static int hammingDistance64(long hash1, long hash2) {
        long i = hash1 ^ hash2;
        i = i - ((i >>> 1) & 0x5555555555555555L);
        i = (i & 0x3333333333333333L) + ((i >>> 2) & 0x3333333333333333L);
        i = (i + (i >>> 4)) & 0x0f0f0f0f0f0f0f0fL;
        i = i + (i >>> 8);
        i = i + (i >>> 16);
        i = i + (i >>> 32);
        return (int) i & 0x7f;
    }
    public int hammingDistance32(int hash1, int hash2) {
        int i = hash1 ^ hash2;
        i = i - ((i >>> 1) & 0x55555555);
        i = (i & 0x33333333) + ((i >>> 2) & 0x33333333);
        i = (i + (i >>> 4)) & 0x0f0f0f0f;
        i = i + (i >>> 8);
        i = i + (i >>> 16);
        return i & 0x3f;
    }

    public static long simhash64(Map<String, Double> tfMap, Map<String, Integer> pepCountMap, int totalPepNum) {
        int bitLen = 64;
        double[] bits = new double[bitLen];
        for (String tag : tfMap.keySet()) {
            double weight = tfMap.get(tag) * Math.log10(totalPepNum/pepCountMap.get(tag)) * 0 +1 ;
            long v = hash64(tag);
            for (int i = bitLen; i >= 1; --i) {
                if (((v >> (bitLen - i)) & 1) == 1)
                    bits[i - 1] += weight;
                else
                    bits[i - 1] -= weight;
            }
        }
        long hash = 0x0000000000000000;
        long one = 0x0000000000000001;
        for (int i = bitLen; i >= 1; --i) {
            if (bits[i - 1] > 0) {
                hash |= one;
            }
            one = one << 1;
        }
        return hash;
    }

    public static int simhash32(Map<String, Double> tfMap, Map<String, Integer> pepCountMap, int totalPepNum) {
        int bitLen = 32;
        double[] bits = new double[bitLen];
//        List<String> tokens = wordSeg.tokens(doc);
        for (String tag : tfMap.keySet()) {
            double weight = tfMap.get(tag) * Math.log10(totalPepNum/pepCountMap.get(tag));
            int v = hash32(tag);
            for (int i = bitLen; i >= 1; --i) {
                if (((v >> (bitLen - i)) & 1) == 1)
                    bits[i - 1] += weight;
                else
                    bits[i - 1] -= weight;
            }
        }
        int hash = 0x00000000;
        int one = 0x00000001;
        for (int i = bitLen; i >= 1; --i) {
            if (bits[i - 1] > 1) {
                hash |= one;
            }
            one = one << 1;
        }
        return hash;
    }
}
