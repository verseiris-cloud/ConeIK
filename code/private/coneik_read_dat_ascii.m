function [t,g,lines,handled,resume]=coneik_read_dat_ascii(fid,strict)
%CONEIK_READ_DAT_ASCII Bounded fast path for ordinary three-column DAT text.
% Four-output callers retain the original rewind-on-failure contract.
% A fifth output permits the caller to retain committed ASCII blocks and
% seek to resume.offset before continuing its unchanged fgetl parser.
% Encoding/other control bytes invalidate the prefix and require a full rewind.
% Raw text/tokens are bounded by two 256 KiB blocks; numeric outputs are O(N).
    t=zeros(0,1);g=zeros(0,2);lines=zeros(0,1);handled=false;
    resume=struct('available',false,'offset',0,'line_count',0,'first',true);
    block_size=256*1024;max_line=64*1024;
    carry=zeros(0,1,'uint8');line_count=0;first=true;
    data_chunks={};line_chunks={};
    while true
        bytes=fread(fid,block_size,'*uint8');at_end=numel(bytes)<block_size;
        % ASCII whitespace 9:13 is accepted by the legacy strtrim/regexp.
        % VT/FF delimit fields, not records; CR/LF splitting stays below.
        % Let fgetl interpret encodings/BOM and all other control bytes.
        if any(bytes>126 | (bytes<32 & (bytes<9 | bytes>13)))
            resume.available=false;return;
        end
        block=[carry;bytes]; %#ok<AGROW>
        if at_end
            last=numel(block);
        else
            endings=find(block==10 | block==13);
            % A terminal CR may be the first half of CRLF in the next block.
            if ~isempty(endings)&&endings(end)==numel(block)&&block(end)==13
                endings(end)=[];
            end
            if isempty(endings),last=0;else,last=endings(end);end
        end
        carry=block(last+1:end);
        % Very long records use the original streaming parser.
        if numel(carry)>max_line,break;end
        if last>0
            raw=regexp(char(block(1:last)'),'\r\n|\n|\r','split')';
            if ~at_end,raw(end)=[];end % Remove the split's synthetic final row.
            if any(cellfun('length',raw)>max_line),break;end
            rowids=line_count+(1:numel(raw))';line_count=line_count+numel(raw);
            raw=strtrim(raw);nonempty=~cellfun('isempty',raw);
            raw=raw(nonempty);rowids=rowids(nonempty);
            if ~isempty(raw)
                tokens=regexp(raw,'[\s,]+','split');sizes=cellfun('length',tokens);
                % Unusual column counts are inexpensive to delegate intact.
                if any(sizes~=3),break;end
                flat=[tokens{:}];values=reshape(str2double(flat),3,[])';
                % A complex token can change storage type for unrelated rows.
                if ~isreal(values),break;end
                if strict
                    if any(~cellfun('isempty',regexp(raw,'(^,)|(,$)|(,\s*,)','once'))),break;end
                end
                numeric=~isnan(values)|reshape(strcmpi(flat,'nan'),3,[])';
                if first && all(~numeric(1,:))
                    values(1,:)=[];rowids(1)=[];
                end
                % Nonfinite/invalid records keep the original mode-specific
                % skip/error behavior, including ordering relative to time checks.
                if any(~isfinite(values),'all'),break;end
                first=false;
                if ~isempty(values)
                    data_chunks{end+1,1}=values;line_chunks{end+1,1}=rowids; %#ok<AGROW>
                end
            end
        end
        % Commit only after every complete record in this block is accepted.
        % Carry starts at this byte position, including a held terminal CR.
        offset=ftell(fid)-numel(carry);
        resume=struct('available',offset>=0,'offset',offset,...
            'line_count',line_count,'first',first);
        if at_end,handled=true;break;end
    end
    if ~handled && (nargout<5 || ~resume.available),return;end
    if ~isempty(data_chunks)
        data=vertcat(data_chunks{:});t=data(:,1);g=data(:,2:3);lines=vertcat(line_chunks{:});
    end
end
