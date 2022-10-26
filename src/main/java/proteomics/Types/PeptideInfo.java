/*
 * Copyright 2016-2019 The Hong Kong University of Science and Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package proteomics.Types;
import java.util.HashSet;
import java.util.Set;

public class PeptideInfo {
    public String seq;
    public boolean isTarget;
    public Set<String> protIdSet = new HashSet<>();
    public char leftFlank;
    public char rightFlank;

    public PeptideInfo(String seq, boolean isTarget, char leftFlank, char rightFlank) {
        this.isTarget = isTarget;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
        this.seq = seq;
    }

    public PeptideInfo clone() throws CloneNotSupportedException {
        super.clone();
        PeptideInfo other = new PeptideInfo(seq, isTarget, leftFlank, rightFlank);
        return other;
    }

}
