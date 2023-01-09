set nocompatible              " be iMproved, required
filetype off                  " required

" set the runtime path to include Vundle and initialize
set rtp+=~/.vim/bundle/Vundle.vim
call vundle#begin()
" alternatively, pass a path where Vundle should install plugins
"call vundle#begin('~/some/path/here')

" let Vundle manage Vundle, required
Plugin 'junegunn/fzf'
Plugin 'VundleVim/Vundle.vim'

Plugin 'scrooloose/nerdtree'
Plugin 'klen/python-mode'
" Plugin 'davidhalter/jedi-vim'

Plugin 'vim-airline/vim-airline'
Plugin 'vim-airline/vim-airline-themes'

"Plugin 'effi/vim-OpenFoam-syntax'
Plugin 'altercation/vim-colors-solarized'
Plugin 'iCyMind/NeoSolarized'
"Plugin 'vim-scripts/xoria256.vim'

Plugin 'majutsushi/tagbar'
"Plugin 'petRUShka/vim-opencl'
Plugin 'w0rp/ale'
Plugin 'digitaltoad/vim-pug'

Plugin 'posva/vim-vue'

Plugin 'lervag/vimtex'
"Plugin 'vim-latex/vim-latex'

Plugin 'Konfekt/FastFold'
Plugin 'matze/vim-tex-fold'
Plugin 'mileszs/ack.vim'

" Track the engine.
Plugin 'SirVer/ultisnips'

" Snippets are separated from the engine. Add this if you want them:
Plugin 'honza/vim-snippets'



Plugin 'pangloss/vim-javascript'
Plugin 'mxw/vim-jsx'
Plugin 'voldikss/vim-translator'
Plugin 'voldikss/vim-floaterm'

Plugin 'leafgarland/typescript-vim'
Plugin 'iamcco/markdown-preview.nvim'
Plugin 'echuraev/translate-shell.vim', { 'do': 'wget -O ~/.vim/trans git.io/trans && chmod +x ~/.vim/trans' }
Plugin 'ruanyl/vim-sort-imports'


" All of your Plugins must be added before the following line
call vundle#end()            " required
filetype plugin indent on    " required
" To ignore plugin indent changes, instead use:
colorschem NeoSolarized
set background=dark
let g:solarized_termcolors=256
"colorscheme solarized
set termguicolors

set t_Co=256
syntax on

set number relativenumber

augroup numbertoggle
  autocmd!
  autocmd BufEnter,FocusGained,InsertLeave * set relativenumber
  autocmd BufLeave,FocusLost,InsertEnter   * set norelativenumber
augroup END

set tabstop=4
set shiftwidth=4
set expandtab
filetype plugin on
set grepprg=grep\ -nH\ $*
filetype indent on
let g:tex_flavor='latex'
let g:Tex_DefaultTargetFormat='pdf'
set linebreak
set showtabline=2

set keymap=russian-jcukenwin


set iminsert=0
set imsearch=0

map <C-n> :NERDTreeToggle<CR>
let NERDTreeIgnore = ['\.pyc$','\.log$','\.aux$','\.blg$']
nmap <F8> :TagbarToggle<CR>


au BufEnter *.tex set keymap=russian-jcukenwin
au BufEnter *.tex setlocal spell spelllang=ru_ru,en
au BufEnter *.tex set tw=80
let g:pymode_rope_lookup_project = 0

nnoremap <buffer> <F9> :exec '!python' shellescape(@%,1)<cr>
nnoremap <F4> :w  <bar> exec '!gfortran '.shellescape('%').' -o '.shellescape('%:r').' && ./'.shellescape('%:r')<CR>
inoremap <A-i> \item

let g:pymode_rope = 0
let g:pymode_rope_completion = 0
let g:pymode_rope_complete_on_dot = 0

" документация
 let g:pymode_doc = 0
 let g:pymode_doc_key = 'K'
" " проверка кода
 let g:pymode_lint = 0
 let g:pymode_lint_checker = "pyflakes,pep8"
 let g:pymode_lint_ignore="E501,W601,C0110"
" " провека кода после сохранения
 let g:pymode_lint_write = 1
"
" " поддержка virtualenv
" let g:pymode_virtualenv = 1
"
" Disable choose first function/method at autocomplete
let g:jedi#popup_select_first = 0
inoremap <C-space> <C-x><C-o>
" text
let g:tex_flavor  = 'latex'
let g:tex_conceal = ''
let g:vimtex_fold_manual = 0
"let g:vimtex_latexmk_continuous = 1
let g:vimtex_compiler_progname = 'nvr'
let g:vimtex_compiler_continuous = 1
"let g:vimtex_latexmk_progname ='nvr'

let g:vimtex_view_method='zathura'
let g:vimtex_quickfix_mode=0
set conceallevel=1
let g:tex_conceal='abdmg'


set nobackup         " no backup files
set nowritebackup    " only in case you don't want a backup file while editing
set noswapfile       " no swap files

" Trigger configuration. Do not use <tab> if you use https://github.com/Valloric/YouCompleteMe.
let g:UltiSnipsExpandTrigger="<tab>"
let g:UltiSnipsJumpForwardTrigger="<c-b>"
let g:UltiSnipsJumpBackwardTrigger="<c-z>"

" If you want :UltiSnipsEdit to split your window.
let g:UltiSnipsEditSplit="vertical"
set linespace=2

" Echo translation in the cmdline
nmap <silent> <Leader>t <Plug>Translate
vmap <silent> <Leader>t <Plug>TranslateV
" Display translation in a window
nmap <silent> <Leader>w <Plug>TranslateW
vmap <silent> <Leader>w <Plug>TranslateWV
" Replace the text with translation
nmap <silent> <Leader>r <Plug>TranslateR
vmap <silent> <Leader>r <Plug>TranslateRV
" Translate the text in clipboard
nmap <silent> <Leader>x <Plug>TranslateX

let g:translator_target_lang = 'en'
let g:translator_source_lang = 'ru'

let g:floaterm_wintype= 'normal'
au FileType javascript setlocal formatprg=prettier
au FileType javascript.jsx setlocal formatprg=prettier
au FileType typescript setlocal formatprg=prettier\ — parser\ typescript

let g:ale_linters = { 'javascript': ['eslint'] }
let g:ale_fixers = { 'javascript': ['eslint'], 'typescript': ['prettier', 'tslint'], 'scss': ['prettier'], 'html': ['prettier'], 'reason': ['refmt'] }
let g:mkdp_filetypes = ['markdown']

g:GuiFont MesloLGS\ NF
let g:vimtex_version_check = 0
