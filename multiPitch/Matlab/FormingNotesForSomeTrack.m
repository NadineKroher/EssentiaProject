function Note = FormingNotesForSomeTrack(track, Group, GroupTrack, F0Location, para)
% Forming notes from the pitch trajectory of some track
%
% Input
%   - track         : track name
%   - Group         : all the notelets
%   - GroupTrack    : track labels of all the notelets
%   - F0Location    : feature vectors of each pitch
%   - para          : parameters
% Output
%   - Note          : notes formed from notelets in this track
%
% Author: Zhiyao Duan
% Created: 5/5/2009
% Last modified: 5/14/2009

mergeNoteGap = para.mergeNoteGap;
minNoteLength = para.minNoteLength;

% initialize notes from groups
groupNum = length(Group);
Note = cell(groupNum, 1);
NoteNum = 0;
for currg = 1:groupNum
    if GroupTrack(currg) == track
        NoteNum = NoteNum + 1;
        Note{NoteNum}(1, :) = F0Location(Group{currg}, 1);  % frame number
        Note{NoteNum}(2, :) = F0Location(Group{currg}, 2);  % frequency
        
        % sort by time
        [Note{NoteNum}(1, :), idx] = sort(Note{NoteNum}(1,:), 'ascend');
        Note{NoteNum}(2, :) = Note{NoteNum}(2, idx);
    end
end
Note = Note(1:NoteNum);

% merge close notes
idx = ones(NoteNum, 1);
for currN = 1:NoteNum
    if isempty(Note{currN}) % since Note{currN} may be changed to 0 in the following
        continue;
    end
    
    for kk = 1:NoteNum
        if kk == currN || isempty(Note{kk})
            continue;
        end
        
        freq = mean(Note{currN}(2,:));          % after merging, these values may change, so re-calculate
        onset = Note{currN}(1,1);
        offset = Note{currN}(1,end);
        
        target_freq = mean(Note{kk}(2,:));
        target_onset = Note{kk}(1,1);
        target_offset = Note{kk}(1,end);
        if abs(target_freq - freq) < 0.5 ...
            && ((target_onset-offset>0 && target_onset-offset<mergeNoteGap)...
                || (target_offset-onset>0 && target_offset-onset<mergeNoteGap))
            if onset < target_onset % connet currN and kk
                fillNote = [offset+1:target_onset-1; ones(1, target_onset-offset-1)*(target_freq + freq)/2];
                Note{currN} = [Note{currN}, fillNote, Note{kk}];
                Note{kk} = [];
                idx(kk) = 0;
            else                    % connet kk and currN
                fillNote = [target_offset+1:onset-1; ones(1, onset-target_offset-1)*(target_freq + freq)/2];
                Note{currN} = [Note{kk}, fillNote, Note{kk}];
                Note{kk} = [];               
                idx(kk) = 0;
            end            
        end
    end
end
% remaining notes are those with idx=1
Note = Note(idx==1);
NoteNum = length(Note);

% remove short notes
idx = ones(NoteNum, 1);
for currN = 1:NoteNum
    if size(Note{currN}, 2) < minNoteLength
        Note{currN} = [];
        idx(currN) = 0;
    end
end
Note = Note(idx==1);